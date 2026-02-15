//! PARSI exhaustive branch-and-bound structural alignment module.
//!
//! Translated from parsi-1.f, parsi-align.f, parsi-score.f, parsi-stack.f, parsi-admin.f.
//!
//! Key notes:
//! - All scoring uses f32 intermediate arithmetic to match Fortran `real`
//! - Distance matrices use f32 sqrt, nint for Fortran-compatible rounding
//! - Expected validation: ~87% on 5-struct, ~78% on 18-struct (search divergence)
//! - ±1 score table diffs from float32 cascade through branch-and-bound

pub mod protein_io;
pub mod scoring;
pub mod domain_tree;
pub mod stack;
pub mod align;

use std::path::Path;
use ndarray::Array2;

use crate::types::ParsiHit;
use crate::store::ProteinStore;
use self::protein_io::{ParsiProtein, parsireadproteindata, treehack, compressca,
                        getdist, hackdist, getdist2, getdist2sum};
use self::scoring::{weights, fillscoretable, selfscore, getupperlower,
                     setcut, flex, setminseglen, setngap, initlexdonestart};
use self::domain_tree::{setldom, init_searchspace, init_ci_ni};
use self::align::align;

// ---------------------------------------------------------------------------
// Constants (from parsizes.for)
// ---------------------------------------------------------------------------

pub const MAXRES1: usize = 4000;
pub const MAXRES2: usize = 4000;
pub const MAXRES0: usize = 400;
pub const EXDIM: usize = 1_000_000;
pub const MAXSTACK: usize = 10_000;
pub const MAXSEG: usize = 200;
pub const MAXDOM: usize = 400;
pub const INFINIT: i32 = 10_000_000;
pub const NUL: i32 = -99;
pub const BOXDIM: usize = 501;
pub const BL: usize = 10;
pub const STARTSIZE: usize = 3;
pub const LENE: usize = 6;
pub const LENH: usize = 8;

// ---------------------------------------------------------------------------
// Cached scoring tables
// ---------------------------------------------------------------------------

use std::sync::OnceLock;

struct ParsiTables {
    weight: Vec<i32>,      // 1001 entries
    scoretable: Vec<i32>,  // 161 * 401 flat
}

static PARSI_TABLES: OnceLock<ParsiTables> = OnceLock::new();

fn get_tables() -> &'static ParsiTables {
    PARSI_TABLES.get_or_init(|| {
        let weight = weights();
        let scoretable = fillscoretable(&weight);
        ParsiTables { weight, scoretable }
    })
}

/// Cached cd1 data for PARSI.
pub struct ParsiCd1Cache {
    pub code: String,
    pub prot: ParsiProtein,
    pub dist: Vec<i16>,       // nres * nres flat
    pub ss: Vec<i32>,         // nseg * nseg flat
    pub upper: Vec<i32>,      // nseg * nseg flat
    pub lower: Vec<i32>,      // nseg * nseg flat
    pub lfix: Vec<bool>,      // nseg * (ndom+1) flat
    pub lfix1: Vec<bool>,     // nseg * (ndom+1) flat
    pub cut: Vec<i32>,        // ndom+1
    pub minseglen: Vec<i32>,  // nseg
    pub ngap: Vec<i32>,       // nseg
    pub dist1sum: Vec<i32>,   // nres * (nres+1) flat
    pub ldom: Vec<bool>,      // ndom+1
}

/// Process a single protein pair through PARSI.
pub fn dowork_parsi(
    cd1: &str,
    cd2: &str,
    dat_dir1: &str,
    dat_dir2: Option<&str>,
    lfirstonly: bool,
    cd1_cache: &mut Option<ParsiCd1Cache>,
) -> Vec<ParsiHit> {
    let tables = get_tables();
    let dat_dir2 = dat_dir2.unwrap_or(dat_dir1);

    let mut output_hits = Vec::new();

    // Phase 1: Load protein 1 (cached if same cd1)
    let need_new = match cd1_cache {
        Some(ref c) => c.code != cd1,
        None => true,
    };

    if need_new {
        let path1 = format!("{}/{}.dat", dat_dir1, cd1);
        let mut p1 = match parsireadproteindata(&path1) {
            Some(p) => p,
            None => return output_hits,
        };

        if p1.nseg <= 2 || p1.ndom < p1.nseg {
            return output_hits;
        }
        if p1.na * p1.na + p1.nb * p1.nb == 0 {
            return output_hits;
        }

        treehack(&mut p1);
        compressca(&mut p1);
        if p1.nres == 0 {
            return output_hits;
        }

        let dist1 = getdist(&p1);
        let ss = selfscore(p1.nseg, &p1.segmentrange, &dist1, p1.nres, &tables.weight);
        let (lower, upper) = getupperlower(&dist1, p1.nseg, &p1.segmentrange, p1.nres);
        let mut dist1_mut = dist1;
        hackdist(&mut dist1_mut, p1.nres);
        let (lfix, lfix1) = flex(&ss, p1.nseg, p1.ndom, &p1.domns, &p1.domseglist);
        let cut = setcut(p1.ndom, &p1.domns, &p1.domseglist, &p1.segmentrange);
        let dist1sum = getdist2sum(&p1.ca, p1.nres);
        let ldom = setldom(p1.ndom, &p1.domns);
        let minseglen = setminseglen(&p1.secstr, p1.nseg);
        let ngap = setngap(&p1.segmentrange, &minseglen, p1.nseg);

        *cd1_cache = Some(ParsiCd1Cache {
            code: cd1.to_string(),
            prot: p1,
            dist: dist1_mut,
            ss, upper, lower,
            lfix, lfix1, cut,
            minseglen, ngap,
            dist1sum, ldom,
        });
    }

    let cache = cd1_cache.as_ref().unwrap();

    // Phase 2: Load protein 2
    let path2 = format!("{}/{}.dat", dat_dir2, cd2);
    let p2 = match parsireadproteindata(&path2) {
        Some(p) => p,
        None => return output_hits,
    };

    if p2.nres == 0 {
        return output_hits;
    }
    if cache.prot.na * p2.na + cache.prot.nb * p2.nb == 0 {
        return output_hits;
    }

    // Build segment mapping for protein 2
    let mut segment2 = vec![0i32; p2.nres];
    for iseg in 0..p2.nseg {
        let s = p2.segmentrange[iseg * 2] as usize;
        let e = p2.segmentrange[iseg * 2 + 1] as usize;
        for ires in (s - 1)..e {
            if ires < p2.nres {
                segment2[ires] = iseg as i32 + 1;
            }
        }
    }

    let dist2sum = getdist2sum(&p2.ca, p2.nres);
    let dist2 = getdist2(&p2.ca, p2.nres);

    // Phase 3: Initialize search space and run alignment
    let nseg = cache.prot.nseg;
    let (trans, mi) = init_searchspace(nseg, p2.nres, &cache.prot.segmentrange,
                                        &cache.minseglen, &cache.ngap, BL);
    let ci0 = init_ci_ni(&mi, nseg);

    let string = format!("{:5}{:5}", cd1, cd2);

    let mut output_lines = Vec::new();
    align(
        &cache.prot, p2.nres, &p2.secstr, p2.nseg, &segment2,
        &cache.dist, cache.prot.nres, &dist2, p2.nres, &dist2sum, &cache.dist1sum,
        &cache.ss, &cache.upper, &cache.lower,
        &cache.prot.segmentrange, &cache.prot.segmentrange0,
        cache.prot.ndom, &cache.prot.node_child, &cache.prot.domns, &cache.prot.domseglist,
        &cache.ldom, &cache.lfix, &cache.lfix1, &cache.cut,
        &cache.minseglen, &cache.ngap,
        &cache.prot.checkrange, &cache.prot.checkx,
        &mi, &ci0, &trans, &tables.scoretable, &tables.weight,
        true, true, lfirstonly, &string, &mut output_lines,
    );

    // Parse output lines into ParsiHit
    for line in &output_lines {
        if let Some(hit) = parse_refine_line(line, cd1, cd2) {
            output_hits.push(hit);
        }
    }

    output_hits
}

/// Parse a refine output line into a ParsiHit.
fn parse_refine_line(line: &str, cd1: &str, cd2: &str) -> Option<ParsiHit> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 4 || !parts[0].starts_with("refine") {
        return None;
    }
    let idom: usize = parts[1].parse().ok()?;
    let score: i32 = parts[2].parse().ok()?;
    let nseg: usize = parts[3].parse().ok()?;
    if parts.len() < 4 + nseg * 4 {
        return None;
    }

    let mut ranges_cd1 = Vec::new();
    let mut ranges_cd2 = Vec::new();
    for i in 0..nseg {
        let a1: i32 = parts[4 + i * 2].parse().ok()?;
        let a2: i32 = parts[4 + i * 2 + 1].parse().ok()?;
        ranges_cd1.push((a1, a2));
    }
    let offset = 4 + nseg * 2;
    for i in 0..nseg {
        let b1: i32 = parts[offset + i * 2].parse().ok()?;
        let b2: i32 = parts[offset + i * 2 + 1].parse().ok()?;
        ranges_cd2.push((b1, b2));
    }

    Some(ParsiHit {
        cd1: cd1.to_string(),
        cd2: cd2.to_string(),
        idom,
        score,
        ranges_cd1,
        ranges_cd2,
    })
}

/// Run PARSI on all pairs from structure lists.
pub fn run_parsi(
    structures: &[String],
    dat_dir: &str,
    dat_dir2: Option<&str>,
) -> Vec<ParsiHit> {
    let _tables = get_tables(); // Initialize tables
    let mut all_hits = Vec::new();
    let mut cd1_cache: Option<ParsiCd1Cache> = None;

    for cd1 in structures {
        for cd2 in structures {
            let hits = dowork_parsi(cd1, cd2, dat_dir, dat_dir2, true, &mut cd1_cache);
            all_hits.extend(hits);
        }
    }

    all_hits
}

/// Run PARSI using ProteinStore (for pipeline integration).
pub fn run_parsi_store(
    structures: &[String],
    store: &ProteinStore,
) -> Vec<ParsiHit> {
    // For now, delegate to file-based version using store's dat_dir
    let dat_dir = store.dat_dir();
    run_parsi(structures, dat_dir, None)
}
