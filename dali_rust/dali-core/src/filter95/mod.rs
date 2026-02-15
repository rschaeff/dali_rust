//! FILTER95 redundancy filtering module.
//!
//! Takes PARSI refine output and:
//! 1. Reads domain tree info (second block, residue-based)
//! 2. Filters by domain type (keep only '*' and '+')
//! 3. Applies FITZ superposition refinement
//! 4. Scores with distance-based scorefun95
//! 5. Filters by Z-score threshold
//!
//! Parameters (from serialcompare.f):
//!   zcut1 = 1.0, fitzrcut = 4.0, fitzmaxiter = 3

use std::collections::HashMap;
use std::sync::Arc;

use ndarray::Array2;

use crate::types::{NodeType, Protein, ParsiHit, DistMatScale100};
use crate::store::ProteinStore;
use crate::numerics::rounding::nint;
use crate::numerics::kabsch::{u3b, transrotate};
use crate::numerics::nw::{filltable_maxsim, nw_maxsim};
use crate::numerics::fitz::FitzResult;

/// FILTER95 output entry.
#[derive(Debug, Clone)]
pub struct Filter95Entry {
    pub cd1: String,
    pub cd2: String,
    pub idom: usize,
    pub score: i32,
    pub zscore: f64,
    pub nseg: usize,
    pub ranges: Vec<i32>, // cd1 ranges then cd2 ranges, nseg*4 values
}

// ---------------------------------------------------------------------------
// Domain info extraction
// ---------------------------------------------------------------------------

/// Domain info cached for a cd1 protein.
struct Cd1DomainInfo {
    lkeep: HashMap<usize, bool>,
    minscore: HashMap<usize, (i32, i32)>,
    segmentrange: HashMap<usize, Vec<(u32, u32)>>,
}

/// Extract domain info from protein's domain tree (second block, residue-based).
fn extract_domain_info(prot: &Protein) -> Cd1DomainInfo {
    let mut lkeep = HashMap::new();
    let mut minscore = HashMap::new();
    let mut segmentrange = HashMap::new();

    for node in &prot.domain_tree {
        let keep = matches!(node.node_type, NodeType::Root | NodeType::Split);
        lkeep.insert(node.index, keep);

        // Compute size from segments (total residue count)
        let size: usize = node
            .segments
            .iter()
            .map(|&(s, e)| (e as usize).saturating_sub(s as usize) + 1)
            .sum();
        minscore.insert(node.index, compute_minscore(size));
        segmentrange.insert(node.index, node.segments.clone());
    }

    Cd1DomainInfo {
        lkeep,
        minscore,
        segmentrange,
    }
}

/// Compute Z-score baseline mean and sigma for a domain size.
///
/// Returns (nint(10000*mean), nint(10000*sigma)).
fn compute_minscore(node_size: usize) -> (i32, i32) {
    let x = (node_size as f64).min(400.0);
    let mean = 7.9494 + 0.70852 * x + 2.5895e-4 * x * x - 1.9156e-6 * x * x * x;
    let sigma = 1.0_f64.max(0.50 * mean);
    (nint(10000.0 * mean), nint(10000.0 * sigma))
}

// ---------------------------------------------------------------------------
// Superposition
// ---------------------------------------------------------------------------

/// Compute rotation/translation from alignment and apply to x.
fn getut_and_transform(x: &mut Array2<f64>, y: &Array2<f64>, ali: &[i32]) {
    let nx = x.ncols();
    let pairs: Vec<(usize, usize)> = (0..nx)
        .filter(|&i| ali[i] >= 0)
        .map(|i| (i, ali[i] as usize))
        .collect();
    let lali = pairs.len();

    if lali == 0 {
        return;
    }

    let mut ux = Array2::zeros((3, lali));
    let mut uy = Array2::zeros((3, lali));
    let w = vec![1.0; lali];

    for (k, &(i, j)) in pairs.iter().enumerate() {
        for dim in 0..3 {
            ux[[dim, k]] = x[[dim, i]];
            uy[[dim, k]] = y[[dim, j]];
        }
    }

    let result = u3b(&w, &ux, &uy, lali, 1);
    transrotate(x, &result.u, &result.t);
}

/// Iterative superposition refinement with initial alignment.
///
/// Unlike numerics::fitz which starts with ali=all -1, this takes
/// an initial alignment and checks convergence against it on the
/// first iteration. Matches Fortran fitz() which takes ali as in-out.
fn fitz95(
    x: &mut Array2<f64>,
    y: &Array2<f64>,
    ali_init: &[i32],
    rcut: f64,
    maxiter: usize,
) -> FitzResult {
    let nx = x.ncols();
    let ny = y.ncols();

    if nx < 1 || ny < 1 {
        return FitzResult {
            ali: vec![-1; nx],
            rms: 0.0,
            lali: 0,
            niter: 0,
        };
    }

    let mut ali = ali_init.to_vec();
    let mut rms = 0.0;
    let mut lali = 0usize;
    let mut niter = 0usize;
    let mut frozen = false;

    while niter < maxiter && !frozen {
        niter += 1;

        let table0 = filltable_maxsim(x, y, rcut);
        let ali_old = ali.clone();
        let (_score, new_ali) = nw_maxsim(nx, ny, &table0);
        ali = new_ali;

        frozen = ali == ali_old;

        // Superpose
        let pairs: Vec<(usize, usize)> = (0..nx)
            .filter(|&i| ali[i] >= 0)
            .map(|i| (i, ali[i] as usize))
            .collect();
        lali = pairs.len();

        if lali > 0 {
            let mut ux = Array2::zeros((3, lali));
            let mut uy = Array2::zeros((3, lali));
            let w = vec![1.0; lali];
            for (k, &(i, j)) in pairs.iter().enumerate() {
                for dim in 0..3 {
                    ux[[dim, k]] = x[[dim, i]];
                    uy[[dim, k]] = y[[dim, j]];
                }
            }
            let result = u3b(&w, &ux, &uy, lali, 1);
            rms = (result.rms / lali.max(1) as f64).sqrt();
            transrotate(x, &result.u, &result.t);
        }
    }

    FitzResult {
        ali,
        rms,
        lali,
        niter,
    }
}

// ---------------------------------------------------------------------------
// Scoring
// ---------------------------------------------------------------------------

/// Distance-based scoring function for FILTER95.
///
/// r1, r2 are int16 distance values (distance * 100).
#[inline]
fn scorefun95(r1: i16, r2: i16) -> f64 {
    let r = (r1 as f64 + r2 as f64) / 200.0;
    let s = (r1 as f64 - r2 as f64) / 100.0;
    if r > 0.01 {
        (0.20 - s.abs() / r) * (-r * r / 400.0).exp()
    } else {
        0.20
    }
}

/// Compute total score from alignment and distance matrices.
///
/// ali is nres1-length, ali[i] > 0 means residue i (0-based) aligns to
/// cd2 residue ali[i] (1-based).
fn gettotscore95(
    ali: &[i32],
    d1: &Array2<i16>,
    nres1: usize,
    d2: &Array2<i16>,
    nres2: usize,
) -> f64 {
    let aligned: Vec<usize> = (0..nres1.min(ali.len()))
        .filter(|&i| ali[i] > 0)
        .collect();

    let mut totscore = 0.0;
    for &k in &aligned {
        for &l in &aligned {
            let q = (ali[k].abs() - 1) as usize;
            let r_idx = (ali[l].abs() - 1) as usize;
            if q < nres2 && r_idx < nres2 {
                totscore += scorefun95(d1[[k, l]], d2[[q, r_idx]]);
            }
        }
    }
    totscore
}

// ---------------------------------------------------------------------------
// Parse PARSI refine lines
// ---------------------------------------------------------------------------

/// Parse a PARSI refine line into a ParsiHit.
///
/// Format: `refine<cd1><cd2> idom score nseg a1 a2 ... b1 b2 ...`
pub fn parse_refine_line(line: &str) -> Option<ParsiHit> {
    let line = line.trim();
    let rpos = line.find("refine")?;

    if rpos + 16 > line.len() {
        return None;
    }

    let cd1 = line[rpos + 6..rpos + 11].to_string();
    let cd2 = line[rpos + 11..rpos + 16].to_string();
    let rest = &line[rpos + 16..];

    let parts: Vec<&str> = rest.split_whitespace().collect();
    if parts.len() < 3 {
        return None;
    }

    let idom: usize = parts[0].parse().ok()?;
    let score: i32 = parts[1].parse().ok()?;
    let nseg: usize = parts[2].parse().ok()?;

    let vals: Vec<i32> = parts[3..]
        .iter()
        .filter_map(|s| s.parse().ok())
        .collect();

    if vals.len() < nseg * 4 {
        return None;
    }

    let mut ranges_cd1 = Vec::with_capacity(nseg);
    let mut ranges_cd2 = Vec::with_capacity(nseg);
    for i in 0..nseg {
        ranges_cd1.push((vals[i * 2], vals[i * 2 + 1]));
        ranges_cd2.push((vals[nseg * 2 + i * 2], vals[nseg * 2 + i * 2 + 1]));
    }

    Some(ParsiHit {
        cd1,
        cd2,
        idom,
        score,
        ranges_cd1,
        ranges_cd2,
    })
}

// ---------------------------------------------------------------------------
// Main FILTER95
// ---------------------------------------------------------------------------

/// Run FILTER95 on ParsiHit results.
///
/// Mirrors comparemodules.f dowork_filter95():
/// - The fast-path Z-score bypass is COMMENTED OUT in Fortran, so ALL
///   alignments go through FITZ refinement.
/// - oldidom is a LOCAL variable initialized to 0 on every call, meaning
///   domain data (resix/xiser/xca) is ALWAYS rebuilt from scratch.
///
/// Parameters:
///   hits: PARSI refine output (ordered by cd1, then cd2)
///   store: protein data store
///   zcut1: Z-score cutoff (default 1.0)
pub fn run_filter95(
    hits: &[ParsiHit],
    store: &ProteinStore,
    zcut1: Option<f64>,
) -> Vec<Filter95Entry> {
    let zcut = zcut1.unwrap_or(1.0);
    let fitzrcut = 4.0;
    let fitzmaxiter = 3;

    let mut results = Vec::new();

    let mut old_cd1 = String::new();
    let mut old_cd2 = String::new();

    // Cached cd1 data
    let mut dinfo: Option<Cd1DomainInfo> = None;
    let mut prot1: Option<Arc<Protein>> = None;
    let mut d1: Option<Arc<DistMatScale100>> = None;

    // Cached cd2 data
    let mut prot2: Option<Arc<Protein>> = None;
    let mut d2: Option<Arc<DistMatScale100>> = None;

    for hit in hits {
        // Load cd1 if changed
        if hit.cd1 != old_cd1 {
            let p = store.get_protein(&hit.cd1).unwrap();
            let dm = store.get_dist_scale100(&hit.cd1).unwrap();
            dinfo = Some(extract_domain_info(&p));
            prot1 = Some(p);
            d1 = Some(dm);
            old_cd1 = hit.cd1.clone();
        }

        let di = dinfo.as_ref().unwrap();
        let p1 = prot1.as_ref().unwrap();
        let dm1 = d1.as_ref().unwrap();

        // Skip non-kept domains.
        // Also correctly rejects idom > ndom (Fortran bug: stale lkeep).
        if !di.lkeep.get(&hit.idom).copied().unwrap_or(false) {
            continue;
        }

        // ALWAYS rebuild domain-local data (oldidom=0 in Fortran module version)
        let segs = match di.segmentrange.get(&hit.idom) {
            Some(s) => s,
            None => continue,
        };

        let mut resix: Vec<usize> = Vec::new();
        let mut xiser: HashMap<usize, usize> = HashMap::new();

        for &(seg_start, seg_end) in segs {
            for j in seg_start as usize..=seg_end as usize {
                if j > 0 && j <= p1.nres {
                    resix.push(j);
                    xiser.insert(j, resix.len()); // 1-based
                }
            }
        }

        let nx = resix.len();
        if nx == 0 {
            continue;
        }

        // Build fresh xca from ca (always rebuilt)
        let mut xca = Array2::<f64>::zeros((3, nx));
        for k in 0..nx {
            let res_idx = resix[k] - 1; // 0-based for ca
            for dim in 0..3 {
                xca[[dim, k]] = p1.ca[[dim, res_idx]];
            }
        }

        // Load cd2 if changed
        if hit.cd2 != old_cd2 {
            let p = store.get_protein(&hit.cd2).unwrap();
            let dm = store.get_dist_scale100(&hit.cd2).unwrap();
            prot2 = Some(p);
            d2 = Some(dm);
            old_cd2 = hit.cd2.clone();
        }

        let p2 = prot2.as_ref().unwrap();
        let dm2 = d2.as_ref().unwrap();

        // Map PARSI alignment to fitzali in domain-local coords
        // fitzali: 0-based domain-local index -> 0-based cd2 residue (or -1)
        let nseg = hit.ranges_cd1.len();
        let mut fitzali = vec![-1i32; nx];

        for i in 0..nseg {
            let (cd1_start, cd1_end) = hit.ranges_cd1[i];
            let (cd2_start, _cd2_end) = hit.ranges_cd2[i];

            if cd1_start > 0 {
                for j in cd1_start..=cd1_end {
                    if let Some(&l) = xiser.get(&(j as usize)) {
                        if l > 0 && l <= nx {
                            // 0-based domain-local index -> 0-based cd2 residue
                            fitzali[l - 1] = cd2_start + j - cd1_start - 1;
                        }
                    }
                }
            }
        }

        // getut + transrotate: initial superposition
        getut_and_transform(&mut xca, &p2.ca, &fitzali);

        // fitz95: iterative refinement with initial alignment
        let fitz_result = fitz95(&mut xca, &p2.ca, &fitzali, fitzrcut, fitzmaxiter);

        // Map back to full-protein alignment for scoring
        // tmpali: nres1-length, tmpali[i] > 0 means 1-based cd2 residue
        let mut tmpali = vec![0i32; p1.nres];
        for i in 0..nx {
            if fitz_result.ali[i] >= 0 {
                let full_res = resix[i]; // 1-based
                tmpali[full_res - 1] = fitz_result.ali[i] + 1; // 1-based
            }
        }

        // Score
        let totscore = gettotscore95(&tmpali, &dm1.data, p1.nres, &dm2.data, p2.nres);
        let xscore = nint(10000.0 * totscore);

        let mut final_score = hit.score;
        let mut final_nseg = nseg;
        let final_ranges: Vec<i32>;

        if xscore > hit.score {
            // Fitz improved the score — use fitz alignment
            final_score = xscore;

            let new_nseg = (0..nx).filter(|&i| fitz_result.ali[i] >= 0).count();
            let mut new_ranges = vec![0i32; new_nseg * 4];
            let mut j = 0;
            for i in 0..nx {
                if fitz_result.ali[i] >= 0 {
                    new_ranges[j * 2] = resix[i] as i32;
                    new_ranges[j * 2 + 1] = resix[i] as i32;
                    new_ranges[j * 2 + new_nseg * 2] = fitz_result.ali[i] + 1;
                    new_ranges[j * 2 + new_nseg * 2 + 1] = fitz_result.ali[i] + 1;
                    j += 1;
                }
            }

            final_nseg = new_nseg;
            final_ranges = new_ranges;
        } else {
            // Use original alignment
            let mut ranges = Vec::with_capacity(nseg * 4);
            for i in 0..nseg {
                ranges.push(hit.ranges_cd1[i].0);
                ranges.push(hit.ranges_cd1[i].1);
            }
            for i in 0..nseg {
                ranges.push(hit.ranges_cd2[i].0);
                ranges.push(hit.ranges_cd2[i].1);
            }
            final_ranges = ranges;
        }

        // Z-score
        let (mean_int, sigma_int) = di.minscore.get(&hit.idom).copied().unwrap_or((0, 1));
        let zscore = (final_score - mean_int) as f64 / sigma_int as f64;

        if zscore < zcut {
            continue;
        }

        results.push(Filter95Entry {
            cd1: hit.cd1.clone(),
            cd2: hit.cd2.clone(),
            idom: hit.idom,
            score: final_score,
            zscore,
            nseg: final_nseg,
            ranges: final_ranges,
        });
    }

    results
}
