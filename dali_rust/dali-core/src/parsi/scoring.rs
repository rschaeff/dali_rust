//! Scoring functions for PARSI — weights, score table, segment scores.

use crate::numerics::rounding::nint;
use super::{NUL, INFINIT, EXDIM, MAXRES0, BL, LENE, LENH};

/// Compute Gaussian envelope weights. w[i] = nint(100*exp(-(i/20)^2)), clipped to 0 if <5.
/// Returns Vec of 1001 entries (index 0 unused, 1..1000). Uses f32.
pub fn weights() -> Vec<i32> {
    let mut w = vec![0i32; 1001];
    let enveloperadius = 20.0f32;
    let x = 1.0f32 / (enveloperadius * enveloperadius);
    for i in 1..1001 {
        let ii = i as f32;
        let val = 100.0f32 * (-x / 100.0f32 * ii * ii).exp();
        w[i] = nint(val as f64) as i32;
        if w[i] < 5 {
            w[i] = 0;
        }
    }
    w
}

/// Compute 160x400 score lookup table. Returns flat Vec indexed [b * 401 + a].
/// scoretable[b, a] = nint(100 * weight[a] * (0.20 - abs(x-y)/x))
pub fn fillscoretable(weight: &[i32]) -> Vec<i32> {
    let mut st = vec![0i32; 161 * 401];
    for b in 1..161 {
        let y = if b < 100 {
            b as f32 * 0.1f32
        } else if b < 125 {
            10.0f32 + (b - 100) as f32 * 0.4f32
        } else {
            (b - 125 + 20) as f32
        };
        for a in 1..401 {
            let x = a as f32 / 10.0f32;
            let val = 100.0f32 * weight[a] as f32 * (0.20f32 - (x - y).abs() / x);
            st[b * 401 + a] = nint(val as f64) as i32;
        }
    }
    st
}

/// Compute segment-segment self-scores using Gaussian envelope.
/// Returns flat Vec of nseg*nseg, indexed [iseg * nseg + jseg].
pub fn selfscore(nseg: usize, segmentrange: &[i32], dist: &[i16], nres: usize, weight: &[i32]) -> Vec<i32> {
    let mut ss = vec![0i32; nseg * nseg];
    let dist_ptr = dist.as_ptr();
    let weight_ptr = weight.as_ptr();
    let weight_len = weight.len();
    for iseg in 0..nseg {
        for jseg in 0..=iseg {
            let si = segmentrange[iseg * 2] as usize;
            let ei = segmentrange[iseg * 2 + 1] as usize;
            let sj = segmentrange[jseg * 2] as usize;
            let ej = segmentrange[jseg * 2 + 1] as usize;
            let mut s: i64 = 0;
            for ires in (si - 1)..ei {
                for jres in (sj - 1)..ej {
                    // Safety: ires < ei <= nres, jres < ej <= nres
                    let p = unsafe { *dist_ptr.add(ires * nres + jres) } as usize;
                    if p < weight_len {
                        s += unsafe { *weight_ptr.add(p) } as i64 * 20;
                    }
                }
            }
            let lali = (ei as i64 - si as i64 + 1) + (ej as i64 - sj as i64 + 1);
            if s < 42 * lali * lali {
                s = 0;
            }
            ss[iseg * nseg + jseg] = s as i32;
            ss[jseg * nseg + iseg] = s as i32;
        }
    }
    ss
}

/// Compute 70%/130% distance-sum bounds per segment pair.
/// Returns (lower, upper) flat Vecs of nseg*nseg.
pub fn getupperlower(dist: &[i16], nseg: usize, segmentrange: &[i32], nres: usize) -> (Vec<i32>, Vec<i32>) {
    let mut lower = vec![0i32; nseg * nseg];
    let mut upper = vec![0i32; nseg * nseg];
    let dist_ptr = dist.as_ptr();
    for iseg in 0..nseg {
        for jseg in 0..=iseg {
            let si = segmentrange[iseg * 2] as usize;
            let ei = segmentrange[iseg * 2 + 1] as usize;
            let sj = segmentrange[jseg * 2] as usize;
            let ej = segmentrange[jseg * 2 + 1] as usize;
            let mut s: i64 = 0;
            for ires in (si - 1)..ei {
                for jres in (sj - 1)..ej {
                    // Safety: ires < ei <= nres, jres < ej <= nres
                    s += unsafe { *dist_ptr.add(ires * nres + jres) } as i64;
                }
            }
            lower[iseg * nseg + jseg] = (70 * s) as i32;
            lower[jseg * nseg + iseg] = (70 * s) as i32;
            upper[iseg * nseg + jseg] = (130 * s) as i32;
            upper[jseg * nseg + iseg] = (130 * s) as i32;
        }
    }
    (lower, upper)
}

/// Compute score cutoffs per domain node.
pub fn setcut(ndom: usize, domns: &[i32], domseglist: &[i32], segmentrange: &[i32]) -> Vec<i32> {
    let alloc = domns.len();
    let mut cut = vec![0i32; alloc];
    for idom in 1..=ndom {
        let ns = domns[idom] as usize;
        let mut lali: i32 = 0;
        for i in 0..ns {
            let seg = domseglist[i * alloc + idom] as usize;
            if seg >= 1 {
                lali += segmentrange[(seg - 1) * 2 + 1] - segmentrange[(seg - 1) * 2] + 1;
            }
        }
        lali = lali.min(200);
        let l = lali as f32;
        let x = 9000.0f32 * (0.83259f32 + 0.11186f32 * l
                              + 3.3537e-5f32 * l * l * l
                              + 1.475e-3f32 * l * l
                              - 1.579e-7f32 * l * l * l * l);
        cut[idom] = x.max(0.0) as i32;
    }
    cut
}

/// Determine fixed/free segments per domain. Returns (lfix, lfix1) flat bool Vecs.
/// lfix[iseg * alloc + idom], lfix1[iseg * alloc + idom].
pub fn flex(ss: &[i32], nseg: usize, ndom: usize, domns: &[i32], domseglist: &[i32]) -> (Vec<bool>, Vec<bool>) {
    let alloc = domns.len();
    let mut lfix = vec![false; nseg * alloc];
    let mut lfix1 = vec![false; nseg * alloc];

    // Compute si and sf
    let mut si = vec![0i32; nseg];
    let mut sf = vec![0.0f32; nseg * nseg];
    for iseg in 0..nseg {
        for jseg in 0..nseg {
            if jseg != iseg {
                si[iseg] += ss[iseg * nseg + jseg];
            }
        }
        if si[iseg] == 0 { si[iseg] = 1; }
        for jseg in 0..nseg {
            if jseg != iseg {
                sf[iseg * nseg + jseg] = ss[iseg * nseg + jseg] as f32 / si[iseg] as f32;
            }
        }
    }

    for idom in (1..=ndom).rev() {
        let ns = domns[idom] as usize;
        let mut seglist = Vec::with_capacity(ns);
        for i in 0..ns {
            seglist.push(domseglist[i * alloc + idom]);
        }
        getlfix(idom, &mut lfix, 0.70, ns, &seglist, &sf, nseg, alloc);
        getlfix(idom, &mut lfix1, 0.90, ns, &seglist, &sf, nseg, alloc);
        // If domain size > 2 and no segment fixed, fix first
        if ns > 2 {
            let n = seglist.iter().filter(|&&seg| {
                let idx = (seg as usize - 1) * alloc + idom;
                idx < lfix.len() && lfix[idx]
            }).count();
            if n == 0 && !seglist.is_empty() {
                let idx = (seglist[0] as usize - 1) * alloc + idom;
                if idx < lfix.len() {
                    lfix[idx] = true;
                }
            }
        }
    }

    (lfix, lfix1)
}

fn getlfix(idom: usize, lfix: &mut Vec<bool>, fixcutoff: f32, ns: usize,
           seglist: &[i32], sf: &[f32], nseg: usize, alloc: usize) {
    for is_ in 0..ns {
        let iseg = (seglist[is_] - 1) as usize; // 0-based
        let mut x = 0.0f32;
        for js in 0..ns {
            if js != is_ {
                let jseg = (seglist[js] - 1) as usize;
                x += sf[iseg * nseg + jseg];
            }
        }
        let idx = iseg * alloc + idom;
        if idx < lfix.len() {
            lfix[idx] = x >= fixcutoff;
        }
    }
}

/// Min segment length: 6 for E, 8 for H.
pub fn setminseglen(secstr: &[char], nseg: usize) -> Vec<i32> {
    let mut minseglen = vec![LENE as i32; nseg];
    for iseg in 0..nseg {
        if iseg < secstr.len() && secstr[iseg] == 'H' {
            minseglen[iseg] = LENH as i32;
        }
    }
    minseglen
}

/// Compute N-terminal gap positions per segment.
pub fn setngap(segmentrange: &[i32], minseglen: &[i32], nseg: usize) -> Vec<i32> {
    let mut ngap = vec![0i32; nseg];
    for iseg in 0..nseg {
        let a2 = segmentrange[iseg * 2 + 1];
        let a1 = segmentrange[iseg * 2];
        let minlen = minseglen[iseg];
        let mut lres = a2 - a1 + 1 - minlen;
        if lres < -29 {
            lres = -29;
        }
        ngap[iseg] = (lres + BL as i32 - 1) / BL as i32;
    }
    ngap
}

/// Initialize lexdone and start arrays.
pub fn initlexdonestart(nseg: usize, ss: &[i32]) -> (Vec<bool>, Vec<i32>) {
    let mut lexdone = vec![false; nseg * nseg];
    let mut start = vec![-INFINIT; nseg * nseg];
    for iseg in 0..nseg {
        for jseg in 0..nseg {
            if ss[iseg * nseg + jseg] <= 0 {
                lexdone[iseg * nseg + jseg] = true;
            }
        }
    }
    (lexdone, start)
}

/// Trim destabilizing end residues from score table.
#[inline(always)]
pub fn trimtable(table: &[i32], table_stride: usize, ibeg: &mut i32, iend: &mut i32,
                  seglen: i32, minlen: i32) -> i32 {
    loop {
        if seglen - *ibeg - *iend <= minlen {
            break;
        }
        // Compute row sums
        let mut rowsum = [0i32; 101];
        for i in (1 + *ibeg)..(seglen - *iend + 1) {
            for j in (1 + *ibeg)..(seglen - *iend + 1) {
                rowsum[i as usize] += table[i as usize * table_stride + j as usize];
            }
        }
        let left = rowsum[(1 + *ibeg) as usize];
        let right = rowsum[(seglen - *iend) as usize];
        if left < right {
            if left < 0 {
                *ibeg += 1;
            } else {
                break;
            }
        } else {
            if right < 0 {
                *iend += 1;
            } else {
                break;
            }
        }
    }

    // Compute total score
    let mut totscore = 0i32;
    for i in (1 + *ibeg)..(seglen - *iend + 1) {
        for j in (1 + *ibeg)..(seglen - *iend + 1) {
            totscore += table[i as usize * table_stride + j as usize];
        }
    }
    totscore
}

/// Compute singlet scores for all candidates of a segment.
pub fn singletex(
    iseg: usize, upp: i32, low: i32, a1: i32, a2: i32,
    nir: usize, trans: &[i32], trans_stride: usize,
    dist2sum: &[i32], dist2sum_stride: usize, nres2: usize,
    dist2: &[i16], nres1: usize, dist: &[i16],
    ex: &mut [i32], bl: usize, start: &[i32], start_stride: usize,
    s_beg: &mut [i16], s_beg_stride: usize,
    s_end: &mut [i16], s_end_stride: usize,
    dist1sum: &[i32], dist1sum_stride: usize,
    minlen: i32, scoretable: &[i32],
) {
    let iwhere = start[iseg * start_stride + iseg];
    if iwhere < 0 { return; }

    let seglen = a2 - a1 + 1;
    if seglen > 100 { return; }

    let dist2_ptr = dist2.as_ptr();
    let dist_ptr = dist.as_ptr();
    let scoretable_ptr = scoretable.as_ptr();

    for ir in 0..nir {
        let transires = trans[ir * trans_stride + iseg];

        if transires == NUL {
            let idx = iwhere as usize + ir;
            if idx < ex.len() { ex[idx] = 0; }
        } else {
            let mut x = -INFINIT;
            for l in 0..bl {
                let tl = transires + l as i32;
                let mut ibeg = 0i32;
                let mut iend = 0i32;

                if tl > nres2 as i32 - minlen + 1 {
                    // Beyond C-terminus
                    if transires != NUL && tl <= nres2 as i32 && tl >= -29 {
                        let key = (tl + 29) as usize;
                        if key < s_beg.len() / s_beg_stride {
                            s_beg[key * s_beg_stride + iseg] = ibeg as i16;
                            s_end[key * s_end_stride + iseg] = iend as i16;
                        }
                    }
                    continue;
                }

                if tl < 1 { ibeg = 1 - tl; }
                if tl + seglen - 1 > nres2 as i32 {
                    iend = seglen - (nres2 as i32 - tl + 1);
                }

                let t;
                let ibeg_save;
                let iend_save;

                if ibeg > seglen - minlen {
                    t = -INFINIT;
                    ibeg_save = INFINIT;
                    iend_save = INFINIT;
                } else {
                    // Calculate per-residue score table (stack-allocated)
                    let tbl_size = 101;
                    let mut table = [0i32; 101 * 101];
                    for i in ibeg..(seglen - iend) {
                        for j in ibeg..(seglen - iend) {
                            let p_row = (a1 - 1 + i) as usize;
                            let p_col = (a1 - 1 + j) as usize;
                            let q_row = (tl - 1 + i) as usize;
                            let q_col = (tl - 1 + j) as usize;
                            if q_row < nres2 && q_col < nres2
                                && p_row < nres1 && p_col < nres1
                            {
                                // Safety: bounds validated by if-guard above
                                let d2_val = unsafe { *dist2_ptr.add(q_row * nres2 + q_col) } as usize;
                                let d1_val = unsafe { *dist_ptr.add(p_row * nres1 + p_col) } as usize;
                                if d2_val < 161 && d1_val < 401 {
                                    table[(i + 1) as usize * tbl_size + (j + 1) as usize] =
                                        unsafe { *scoretable_ptr.add(d2_val * 401 + d1_val) };
                                }
                            }
                        }
                    }

                    t = trimtable(&table, tbl_size, &mut ibeg, &mut iend, seglen, minlen);
                    ibeg_save = ibeg;
                    iend_save = iend;
                }

                // Remember end-gaps
                if transires != NUL && tl <= nres2 as i32 && tl >= -29 {
                    let key = (tl + 29) as usize;
                    if key < s_beg.len() / s_beg_stride {
                        s_beg[key * s_beg_stride + iseg] = if ibeg <= seglen { ibeg_save as i16 } else { INFINIT as i16 };
                        s_end[key * s_end_stride + iseg] = if iend <= seglen { iend_save as i16 } else { INFINIT as i16 };
                    }
                }

                x = x.max(t);
            }

            let idx = iwhere as usize + ir;
            if idx < ex.len() { ex[idx] = x; }
        }
    }
}

/// Compute segment-segment doublet score.
#[inline(always)]
pub fn segsegscore(
    iseg: usize, jseg: usize, transires: i32, transjres: i32,
    a1: i32, a2: i32, b1: i32, b2: i32,
    dist: &[i16], nres1: usize, dist2: &[i16], nres2: usize,
    upp1: i32, low1: i32,
    dist2sum: &[i32], dist2sum_stride: usize,
    bl: usize, lseqtl: bool,
    s_beg: &[i16], s_beg_stride: usize,
    s_end: &[i16], s_end_stride: usize,
    dist1sum: &[i32], dist1sum_stride: usize,
    scoretable: &[i32],
) -> i32 {
    if transires == NUL || transjres == NUL { return 0; }

    let mut s = -INFINIT;
    if lseqtl && iseg > jseg && transires < transjres { return s; }
    if lseqtl && iseg < jseg && transires > transjres { return s; }

    let lself = iseg == jseg;
    let j1 = transires - a1;
    let j2 = transjres - b1;

    let dist2_ptr = dist2.as_ptr();
    let dist_ptr = dist.as_ptr();
    let scoretable_ptr = scoretable.as_ptr();

    for jshift in 0..bl {
        if transjres + jshift as i32 > nres2 as i32 { break; }
        let jbeg_key = transjres + jshift as i32;
        if jbeg_key < -29 || jbeg_key > nres2 as i32 { continue; }
        let key_j = (jbeg_key + 29) as usize;
        if key_j >= s_beg.len() / s_beg_stride { continue; }
        let jbeg = s_beg[key_j * s_beg_stride + jseg] as i32;
        let jend = s_end[key_j * s_end_stride + jseg] as i32;
        let l0 = jbeg == 0 && jend == 0;
        let k2 = j2 + jshift as i32;
        let d1_val = b1 + k2 - 1 + jbeg;
        let d2_val = b2 + k2 - jend;
        let e1 = b1 + jbeg - 1;
        let e2 = b2 - jend;

        let jfirst = transjres + jshift as i32 + jbeg;
        let jlast = transjres + jshift as i32 - jend + b2 - b1;
        if jlast > nres2 as i32 { break; }

        for ishift in 0..bl {
            if lself && jshift != ishift { continue; }
            if transires + ishift as i32 > nres2 as i32 { continue; }

            let ibeg_key = transires + ishift as i32;
            if ibeg_key < -29 || ibeg_key > nres2 as i32 { continue; }
            let key_i = (ibeg_key + 29) as usize;
            if key_i >= s_beg.len() / s_beg_stride { continue; }
            let ibeg = s_beg[key_i * s_beg_stride + iseg] as i32;
            let iend = s_end[key_i * s_end_stride + iseg] as i32;
            let l1 = l0 && ibeg == 0 && jbeg == 0;
            let k1 = j1 + ishift as i32;

            let ifirst = transires + ishift as i32 + ibeg;
            let ilast = transires + ishift as i32 - iend + a2 - a1;
            if ilast > nres2 as i32 { continue; }

            // Disallow segment overlaps
            if iseg < jseg && ilast >= jfirst { continue; }
            if jseg < iseg && jlast >= ifirst { continue; }

            // Distance sum filter
            let mut ds2: i64 = 0;
            for c in (a1 + k1 + ibeg)..=(a2 + k1 - iend) {
                if c >= 1 && c <= nres2 as i32 {
                    let ci = (c - 1) as usize;
                    ds2 += (dist2sum[ci * dist2sum_stride + d2_val as usize]
                            - dist2sum[ci * dist2sum_stride + d1_val as usize]) as i64;
                }
            }

            if !l1 {
                let mut ds1: i64 = 0;
                for c in (a1 + ibeg)..=(a2 - iend) {
                    if c >= 1 && c <= nres1 as i32 {
                        let ci = (c - 1) as usize;
                        ds1 += (dist1sum[ci * dist1sum_stride + e2 as usize]
                                - dist1sum[ci * dist1sum_stride + e1 as usize]) as i64;
                    }
                }
                let low = (0.7f32 * ds1 as f32) as i64;
                let upp = (1.3f32 * ds1 as f32) as i64;
                if ds2 <= low || ds2 >= upp { continue; }
            } else {
                if ds2 <= low1 as i64 || ds2 >= upp1 as i64 { continue; }
            }

            // Precise calculation — unsafe unchecked indexing on hot path
            let mut x: i32 = 0;
            for a in (a1 + ibeg)..=(a2 - iend) {
                let p_row = (a - 1) as usize;
                for b in (b1 + jbeg)..=(b2 - jend) {
                    let b_col = (b - 1) as usize;
                    let d2_row = (transires + ishift as i32 - 1 + (a - a1)) as usize;
                    let d2_col = (transjres + jshift as i32 - 1 + (b - b1)) as usize;
                    if d2_row < nres2 && d2_col < nres2
                        && p_row < nres1 && b_col < nres1
                    {
                        // Safety: bounds validated by if-guard above;
                        // dv2 < 161 && dv1 < 401 ensures scoretable index < 161*401
                        let dv2 = unsafe { *dist2_ptr.add(d2_row * nres2 + d2_col) } as usize;
                        let dv1 = unsafe { *dist_ptr.add(p_row * nres1 + b_col) } as usize;
                        if dv2 >= 1 && dv1 >= 1 && dv2 < 161 && dv1 < 401 {
                            x += unsafe { *scoretable_ptr.add(dv2 * 401 + dv1) };
                        }
                    }
                }
            }
            if x > s { s = x; }
        }
    }

    s
}

/// Compute doublet scores for all candidate pairs of two segments.
pub fn doubletex(
    aseg: usize, bseg: usize, nir: usize, njr: usize,
    upp: i32, low: i32, a1: i32, a2: i32, b1: i32, b2: i32,
    trans: &[i32], trans_stride: usize,
    dist2sum: &[i32], dist2sum_stride: usize, nres2: usize,
    dist2: &[i16], nres1: usize, dist: &[i16],
    ex: &mut [i32], bl: usize,
    start: &[i32], start_stride: usize, lseqtl: bool,
    s_beg: &[i16], s_beg_stride: usize,
    s_end: &[i16], s_end_stride: usize,
    dist1sum: &[i32], dist1sum_stride: usize,
    scoretable: &[i32],
) {
    let ijstart = start[aseg * start_stride + bseg];
    if ijstart < 0 { return; }

    for ir in 0..nir {
        let iwhere = ijstart + ir as i32 * njr as i32;
        let transires = trans[ir * trans_stride + aseg];
        if transires == NUL {
            for jr in 0..njr {
                let idx = (iwhere + jr as i32) as usize;
                if idx < ex.len() { ex[idx] = 0; }
            }
        } else {
            for jr in 0..njr {
                let transjres = trans[jr * trans_stride + bseg];
                let x = if transjres == NUL {
                    0
                } else {
                    segsegscore(
                        aseg, bseg, transires, transjres,
                        a1, a2, b1, b2, dist, nres1, dist2, nres2,
                        upp, low, dist2sum, dist2sum_stride,
                        bl, lseqtl,
                        s_beg, s_beg_stride, s_end, s_end_stride,
                        dist1sum, dist1sum_stride, scoretable,
                    )
                };
                let idx = (iwhere + jr as i32) as usize;
                if idx < ex.len() { ex[idx] = x; }
            }
        }
    }
}

/// Compute scores for new segment pairs in domain idom.
pub fn update_ex(
    mi: &[i32], trans: &[i32], trans_stride: usize,
    upper: &[i32], lower: &[i32], nseg: usize,
    segmentrange: &[i32],
    idom: usize, domns: &[i32], domseglist: &[i32], alloc: usize,
    ex: &mut [i32],
    dist: &[i16], nres1: usize, dist2: &[i16], nres2: usize,
    dist2sum: &[i32], dist2sum_stride: usize,
    nix: &mut usize, laststart: &mut i32,
    ixstart: &mut [i32],
    bl: usize, start: &mut [i32],
    lexdone: &mut [bool], lseqtl: bool,
    s_beg: &mut [i16], s_beg_stride: usize,
    s_end: &mut [i16], s_end_stride: usize,
    dist1sum: &[i32], dist1sum_stride: usize,
    minseglen: &[i32], scoretable: &[i32],
) {
    let ns = domns[idom] as usize;
    for is_ in 0..ns {
        let iseg = (domseglist[is_ * alloc + idom] - 1) as usize;
        if !lexdone[iseg * nseg + iseg] {
            setstart1(iseg, iseg, nix, laststart, ixstart, mi, start, nseg, lexdone);
            singletex(
                iseg, upper[iseg * nseg + iseg], lower[iseg * nseg + iseg],
                segmentrange[iseg * 2], segmentrange[iseg * 2 + 1],
                mi[iseg] as usize, trans, trans_stride,
                dist2sum, dist2sum_stride, nres2,
                dist2, nres1, dist, ex, bl, start, nseg,
                s_beg, s_beg_stride, s_end, s_end_stride,
                dist1sum, dist1sum_stride,
                minseglen[iseg], scoretable,
            );
        }

        for js in 0..is_ {
            let jseg = (domseglist[js * alloc + idom] - 1) as usize;
            if !lexdone[iseg * nseg + jseg] {
                let aseg = iseg.min(jseg);
                let bseg = iseg.max(jseg);
                setstart1(aseg, bseg, nix, laststart, ixstart, mi, start, nseg, lexdone);
                doubletex(
                    aseg, bseg, mi[aseg] as usize, mi[bseg] as usize,
                    upper[aseg * nseg + bseg], lower[aseg * nseg + bseg],
                    segmentrange[aseg * 2], segmentrange[aseg * 2 + 1],
                    segmentrange[bseg * 2], segmentrange[bseg * 2 + 1],
                    trans, trans_stride,
                    dist2sum, dist2sum_stride, nres2,
                    dist2, nres1, dist, ex, bl,
                    start, nseg, lseqtl,
                    s_beg, s_beg_stride, s_end, s_end_stride,
                    dist1sum, dist1sum_stride, scoretable,
                );
            }
        }
    }
}

fn setstart1(
    iseg: usize, jseg: usize,
    nix: &mut usize, laststart: &mut i32,
    ixstart: &mut [i32], mi: &[i32],
    start: &mut [i32], nseg: usize,
    lexdone: &mut [bool],
) {
    ixstart[0 * (nseg * nseg) + *nix] = iseg as i32;
    ixstart[1 * (nseg * nseg) + *nix] = jseg as i32;
    start[iseg * nseg + jseg] = *laststart;
    start[jseg * nseg + iseg] = *laststart;
    if iseg == jseg {
        *laststart += mi[iseg];
    } else {
        *laststart += mi[iseg] * mi[jseg];
    }
    if *laststart > EXDIM as i32 {
        *laststart = 0;
    }
    lexdone[iseg * nseg + jseg] = true;
    lexdone[jseg * nseg + iseg] = true;
    *nix += 1;
}

/// Compute segment score contributions (ess) and total estimate.
/// Writes into caller-provided `ess` buffer (nseg*nseg).
#[inline(always)]
pub fn get_ess_into(
    ns: usize, seglist: &[usize], ni: &[i32], ci: &[i32], ci_stride: usize,
    ex: &[i32], start: &[i32], start_stride: usize,
    mi: &[i32], trans: &[i32], trans_stride: usize, lseqtl: bool,
    ess: &mut [i32],
) -> i32 {
    let nseg = mi.len();
    ess[..nseg * nseg].fill(0);
    let mut est: i64 = 0;
    for is_ in 0..ns {
        let iseg = seglist[is_];
        let e = get_estimate(iseg, iseg, ni, ci, ci_stride, ex, start, start_stride,
                              mi, trans, trans_stride, lseqtl);
        ess[iseg * nseg + iseg] = e;
        est += e as i64;
        for js in 0..is_ {
            let jseg = seglist[js];
            let e = get_estimate(iseg, jseg, ni, ci, ci_stride, ex, start, start_stride,
                                  mi, trans, trans_stride, lseqtl);
            let e2 = e + e;
            ess[iseg * nseg + jseg] = e2;
            ess[jseg * nseg + iseg] = e2;
            est += e2 as i64;
        }
    }
    est.min(i32::MAX as i64) as i32
}

/// Convenience wrapper that allocates ess buffer.
#[inline(always)]
pub fn get_ess(
    ns: usize, seglist: &[usize], ni: &[i32], ci: &[i32], ci_stride: usize,
    ex: &[i32], start: &[i32], start_stride: usize,
    mi: &[i32], trans: &[i32], trans_stride: usize, lseqtl: bool,
) -> (Vec<i32>, i32) {
    let nseg = mi.len();
    let mut ess = vec![0i32; nseg * nseg];
    let est = get_ess_into(ns, seglist, ni, ci, ci_stride, ex, start, start_stride,
                            mi, trans, trans_stride, lseqtl, &mut ess);
    (ess, est)
}

/// Get maximum score estimate for a segment pair.
#[inline(always)]
pub fn get_estimate(
    iseg: usize, jseg: usize,
    ni: &[i32], ci: &[i32], ci_stride: usize,
    ex: &[i32], start: &[i32], start_stride: usize,
    mi: &[i32], trans: &[i32], trans_stride: usize, lseqtl: bool,
) -> i32 {
    let aseg = iseg.min(jseg);
    let bseg = iseg.max(jseg);

    if start[iseg * start_stride + jseg] < 0 {
        let mut est = 0;
        if lseqtl && aseg != bseg && ni[aseg] >= 1 && ni[bseg] >= 1 {
            let ares = trans[(ci[(ni[aseg] - 1) as usize * ci_stride + aseg] as usize) * trans_stride + aseg];
            if ares == NUL { return 0; }
            let bres = trans[(ci[(ni[bseg] - 1) as usize * ci_stride + bseg] as usize) * trans_stride + bseg];
            if bres == NUL { return 0; }
            let ares0 = trans[(ci[0 * ci_stride + aseg] as usize) * trans_stride + aseg];
            if ares0 > bres { est = -INFINIT; }
        }
        return est;
    }

    let mut est = -INFINIT;
    let iwhere0 = start[aseg * start_stride + bseg];
    let ex_ptr = ex.as_ptr();
    let ex_len = ex.len();
    if aseg == bseg {
        for i in 0..ni[aseg] as usize {
            let idx = iwhere0 as usize + ci[i * ci_stride + aseg] as usize;
            if idx < ex_len {
                // Safety: idx < ex_len verified above
                let val = unsafe { *ex_ptr.add(idx) };
                if val > est { est = val; }
            }
        }
    } else {
        let jw = mi[bseg];
        for i in 0..ni[aseg] as usize {
            let iwhere = iwhere0 + ci[i * ci_stride + aseg] * jw;
            for j in 0..ni[bseg] as usize {
                let idx = (iwhere + ci[j * ci_stride + bseg]) as usize;
                if idx < ex_len {
                    // Safety: idx < ex_len verified above
                    let val = unsafe { *ex_ptr.add(idx) };
                    if val > est { est = val; }
                }
            }
        }
    }
    est
}
