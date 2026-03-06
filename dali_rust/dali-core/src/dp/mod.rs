//! DP module: Z-score computation from WOLF alignments.
//!
//! Translated from comparemodules.f dp module.

use ndarray::Array2;

use crate::types::{
    AlignmentBlock, DccpEntry, DomainNode,
};
use crate::numerics::scoring::{dpweights, totscore, zscore_func};
use crate::numerics::calc_rmsd;
use crate::store::ProteinStore;

/// Default Z-score cutoff.
pub const ZCUT_DEFAULT: f64 = 2.0;

/// Score a single domain pair.
///
/// Builds alignment array from blocks, filtering by domain active residues,
/// then computes raw score via totscore() and Z-score via zscore_func().
///
/// Translated from comparemodules.f dopair().
///
/// # Arguments
/// * `idom1`, `idom2` - domain indices into filtered_domains lists
/// * `blocks` - alignment blocks (1-based inclusive ranges)
/// * `nres1`, `nres2` - residue counts
/// * `d1`, `d2` - distance matrices (scale 10)
/// * `doms1`, `doms2` - filtered domain node lists
///
/// # Returns
/// (zscore, raw_score)
pub fn dopair(
    idom1: usize,
    idom2: usize,
    blocks: &[AlignmentBlock],
    nres1: usize,
    nres2: usize,
    d1: &Array2<i16>,
    d2: &Array2<i16>,
    doms1: &[&DomainNode],
    doms2: &[&DomainNode],
    wght: &[f64; 101],
) -> (f64, f64) {
    // Mark active residues in each domain (1-based indexing)
    let mut lactive1 = vec![false; nres1 + 1];
    let mut lactive2 = vec![false; nres2 + 1];

    let mut len1 = 0usize;
    for &(seg_start, seg_end) in &doms1[idom1].segments {
        for j in seg_start..=seg_end.min(nres1 as u32) {
            lactive1[j as usize] = true;
            len1 += 1;
        }
    }

    let mut len2 = 0usize;
    for &(seg_start, seg_end) in &doms2[idom2].segments {
        for j in seg_start..=seg_end.min(nres2 as u32) {
            lactive2[j as usize] = true;
            len2 += 1;
        }
    }

    // Build alignment array: ali1[i] = j (1-based) or 0
    let mut ali1 = vec![0i32; nres1];
    let mut lali = 0usize;

    for b in blocks {
        for j in b.l1..=b.r1 {
            let j2 = b.l2 as i32 + j as i32 - b.l1 as i32;
            let j_u = j as usize;
            let j2_u = j2 as usize;
            if j_u >= 1 && j_u <= nres1 && j2_u >= 1 && j2_u <= nres2 {
                if lactive1[j_u] && lactive2[j2_u] {
                    ali1[j_u - 1] = j2 as i32; // 1-based target in 0-based array
                    lali += 1;
                }
            }
        }
    }

    if lali == 0 {
        return (0.0, 0.0);
    }

    let x = totscore(&ali1, nres1, d1, d2, wght);
    let z = zscore_func(len1, len2, x);

    (z, x)
}

/// Run DP scoring on a list of alignments (from WOLF).
///
/// For each alignment, loads proteins from the store, scores all domain pairs,
/// outputs a DccpEntry if any domain pair has Z-score >= zcut.
///
/// Matches the Fortran dowork_dp loop.
///
/// # Arguments
/// * `alignments` - list of (cd1, cd2, blocks) from WOLF
/// * `store` - protein store for loading data
/// * `zcut` - Z-score cutoff (default 2.0)
pub fn run_dp(
    alignments: &[(String, String, Vec<AlignmentBlock>)],
    store: &ProteinStore,
    zcut: f64,
) -> Vec<DccpEntry> {
    let wght = dpweights();
    let mut results = Vec::new();

    for (cd1, cd2, blocks) in alignments {
        // Load proteins and distance matrices
        let prot1 = match store.get_protein(cd1) {
            Ok(p) => p,
            Err(_) => continue,
        };
        let prot2 = match store.get_protein(cd2) {
            Ok(p) => p,
            Err(_) => continue,
        };
        let d1 = match store.get_dist_scale10(cd1) {
            Ok(d) => d,
            Err(_) => continue,
        };
        let d2 = match store.get_dist_scale10(cd2) {
            Ok(d) => d,
            Err(_) => continue,
        };

        let doms1 = prot1.filtered_domains();
        let doms2 = prot2.filtered_domains();

        if doms1.is_empty() || doms2.is_empty() {
            continue;
        }

        // Score all domain pairs
        let mut zmax = 0.0f64;
        let mut x1 = 0.0f64;

        for idom1 in 0..doms1.len() {
            for idom2 in 0..doms2.len() {
                let (z, x) = dopair(
                    idom1, idom2, blocks,
                    prot1.nres, prot2.nres,
                    &d1.data, &d2.data,
                    &doms1, &doms2,
                    &wght,
                );
                if z > zmax {
                    zmax = z;
                }
                if idom1 == 0 && idom2 == 0 {
                    x1 = x;
                }
            }
        }

        // Output if any domain pair passed z cutoff
        if zmax >= zcut {
            let rmsd = calc_rmsd(&prot1.ca, &prot2.ca, blocks);
            results.push(DccpEntry {
                cd1: cd1.clone(),
                cd2: cd2.clone(),
                score: x1,
                zscore: zmax,
                rmsd,
                blocks: blocks.clone(),
            });
        }
    }

    results
}
