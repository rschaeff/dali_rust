//! Domain tree handling for PARSI.

use super::{NUL, BL, MAXRES0, STARTSIZE};

/// Mark domains with >= STARTSIZE segments as 'align'.
/// Returns boolean vec ldom[0..=ndom], with ldom[0]=false.
pub fn setldom(ndom: usize, domns: &[i32]) -> Vec<bool> {
    let mut ldom = vec![false; ndom + 1];
    for idom in 1..=ndom {
        if idom < domns.len() {
            ldom[idom] = domns[idom] >= STARTSIZE as i32;
        }
    }
    ldom
}

/// Initialize search space: enumerate candidate residues for each segment.
///
/// For each segment, candidates are:
/// 1. Residues in protein 2 at bl-residue steps
/// 2. N-terminal gap residues (negative indices down to -29)
/// 3. NUL candidate (unaligned)
///
/// Returns (trans, mi) where:
///   trans: flat [max_cand * nseg], indexed [ir * nseg + iseg]
///   mi: [nseg] — number of candidates per segment
pub fn init_searchspace(
    nseg: usize, nres2: usize, segmentrange: &[i32],
    minseglen: &[i32], ngap: &[i32], bl: usize,
) -> (Vec<i32>, Vec<i32>) {
    let max_cand = MAXRES0;
    let mut trans = vec![NUL; max_cand * nseg];
    let mut mi = vec![0i32; nseg];

    for iseg in 0..nseg {
        let mut ir: usize = 0;

        // Regular candidates: step through protein 2 at bl intervals
        let mut ires = 1i32;
        while ires <= nres2 as i32 {
            if ir >= max_cand { break; }
            trans[ir * nseg + iseg] = ires;
            ir += 1;
            ires += bl as i32;
        }

        // N-terminal gap candidates
        for igap in 1..=ngap[iseg] {
            if ir >= max_cand { break; }
            let mut gapres = 1 - igap * bl as i32;
            if gapres < -29 { gapres = -29; }
            trans[ir * nseg + iseg] = gapres;
            ir += 1;
        }

        // NUL candidate
        if ir < max_cand {
            trans[ir * nseg + iseg] = NUL;
            ir += 1;
        }

        mi[iseg] = ir as i32;
    }

    (trans, mi)
}

/// Initialize candidate index arrays — all candidates active.
/// ci[ir * nseg + iseg] = ir (identity mapping)
/// Returns ci (flat [max_cand * nseg]).
pub fn init_ci_ni(mi: &[i32], nseg: usize) -> Vec<i32> {
    let max_cand = MAXRES0;
    let mut ci = vec![0i32; max_cand * nseg];
    for iseg in 0..nseg {
        for ir in 0..mi[iseg] as usize {
            ci[ir * nseg + iseg] = ir as i32;
        }
    }
    ci
}
