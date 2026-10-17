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
/// Port of parsi-align.f init_searchspace (non-prealigned branch; lpreali is a
/// .false. parameter in comparemodules.f). Per segment, in this order:
/// 1. N-terminal gap starts: `ires=1; do i=ngap,1,-1: ires=ires-i*bl` -- the
///    decrement is CUMULATIVE, so ngap=2 gives -19,-29 and ngap=3 gives -29,-49,-59
///    (starts below -29 can never score; they still occupy candidate slots)
/// 2. Residues in protein 2 at bl-residue steps
/// 3. NUL (lnul is .true.)
///
/// Returns (trans, mi) where:
///   trans: flat [max_cand * nseg], indexed [ir * nseg + iseg]
///   mi: [nseg] -- number of candidates per segment
pub fn init_searchspace(
    nseg: usize, nres2: usize, _segmentrange: &[i32],
    _minseglen: &[i32], ngap: &[i32], bl: usize,
) -> (Vec<i32>, Vec<i32>) {
    let max_cand = MAXRES0;
    let mut trans = vec![NUL; max_cand * nseg];
    let mut mi = vec![0i32; nseg];

    for iseg in 0..nseg {
        let mut ir: usize = 0;
        let mut push = |trans: &mut Vec<i32>, ir: &mut usize, v: i32| {
            if *ir < max_cand {
                trans[*ir * nseg + iseg] = v;
                *ir += 1;
            }
        };

        let mut ires = 1i32;
        let mut i = ngap[iseg];
        while i >= 1 {
            ires -= i * bl as i32;
            push(&mut trans, &mut ir, ires);
            i -= 1;
        }
        let mut r = 1i32;
        while r <= nres2 as i32 {
            push(&mut trans, &mut ir, r);
            r += bl as i32;
        }
        push(&mut trans, &mut ir, NUL);

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
