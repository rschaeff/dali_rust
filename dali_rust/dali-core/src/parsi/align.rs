//! Core alignment module for PARSI — align() orchestrator + helpers.

use std::collections::HashMap;
use super::{NUL, INFINIT, EXDIM, MAXRES0, BL};
use super::scoring::{get_ess, get_estimate, update_ex, initlexdonestart,
                      singletex, doubletex};
use super::stack::{PriorityStack, SearchState, getnextbest, declump};

/// Main PARSI alignment — bottom-up over domain tree.
#[allow(clippy::too_many_arguments)]
pub fn align(
    prot1: &super::protein_io::ParsiProtein,
    prot2_nres: usize, prot2_secstr: &[char], prot2_nseg: usize, prot2_segment: &[i32],
    dist: &[i16], nres1: usize, dist2: &[i16], nres2: usize,
    dist2sum: &[i32], dist1sum: &[i32],
    ss: &[i32], upper: &[i32], lower: &[i32],
    segmentrange: &[i32], segmentrange0: &[i32],
    ndom: usize, node_child: &[i32], domns: &[i32], domseglist: &[i32],
    ldom: &[bool], lfix: &[bool], lfix1: &[bool], cut: &[i32],
    minseglen: &[i32], ngap: &[i32],
    checkrange: &[i32], checkx: &[i32],
    mi: &[i32], ci0: &[i32], trans: &[i32],
    scoretable: &[i32], weight: &[i32],
    lseqtl: bool, lnul: bool, lfirstonly: bool,
    string: &str, outputunit_lines: &mut Vec<String>,
) {
    let nseg = prot1.nseg;
    let alloc = node_child.len() / 2; // ndom_alloc

    // Initialize ex array and tracking
    let (mut lexdone, mut start) = initlexdonestart(nseg, ss);
    let mut ex = vec![0i32; EXDIM];
    let mut nix: usize = 0;
    let mut laststart: i32 = 0;
    let mut ixstart = vec![0i32; 2 * nseg * nseg];

    // Per-residue score trimming arrays
    let s_beg_stride = nseg;
    let s_end_stride = nseg;
    let s_beg_rows = nres2 + 30;
    let mut s_beg = vec![0i16; s_beg_rows * s_beg_stride];
    let mut s_end = vec![0i16; s_beg_rows * s_end_stride];

    // Ali save arrays for child results
    let ali_save_len = 100001;
    let mut ali_save = vec![0i32; ali_save_len];
    let ali_tree_len = ndom + 100;
    let mut ali_start = vec![0i32; ali_tree_len];
    let mut ali_nali = vec![0i32; ali_tree_len];
    let mut ali_laststart: i32 = 0;

    // Ali save for per-residue (findhomolog)
    let mut ali_save1 = vec![0i32; ali_save_len];
    let mut ali_start1 = vec![0i32; ali_tree_len];
    let mut ali_nali1 = vec![0i32; ali_tree_len];
    let mut ali_laststart1: i32 = 0;

    // Working copy of trans for shipup
    let mut trans1 = vec![0i32; trans.len()];
    let mut ni = vec![0i32; nseg];

    // Mutable copies of mi, ci0, trans for compress to modify
    let mut mi_w = mi.to_vec();
    let mut ci0_w = ci0.to_vec();
    let mut trans_w = trans.to_vec();

    // Cut working copy
    let mut cut_w = cut.to_vec();

    let dist2sum_stride = nres2 + 1;
    let dist1sum_stride = nres1 + 1;

    // Process domain tree bottom-up
    for idom in (1..=ndom).rev() {
        let mut icyc = 0;
        if idom >= ldom.len() || !ldom[idom] { continue; }

        let idom1 = node_child[0 * alloc + idom] as usize;
        let idom2 = node_child[1 * alloc + idom] as usize;
        let ns = domns[idom] as usize;
        let mut seglist_1 = Vec::with_capacity(ns);
        for i in 0..ns {
            seglist_1.push(domseglist[i * alloc + idom]);
        }
        // 0-based segment indices
        let seglist: Vec<usize> = seglist_1.iter().map(|&s| (s - 1) as usize).collect();

        // shipup: build work search space from child node results
        shipup(&mut ni, &mut trans1, idom, ldom, domns, domseglist, alloc,
               &mi_w, &ci0_w, &trans_w, idom1, idom2,
               &ali_start, &ali_nali, &ali_save, BL, nres2, lnul);

        // compress: remap ex[] indices after merging child nodes
        if idom1 > 0 {
            compress(idom1, &ni, &trans1, &mut mi_w, &mut ci0_w, &mut trans_w,
                     &mut ex, &mut nix, &mut laststart, &mut ixstart,
                     domns, domseglist, alloc, &mut start, nseg);
        }
        if idom2 > 0 {
            compress(idom2, &ni, &trans1, &mut mi_w, &mut ci0_w, &mut trans_w,
                     &mut ex, &mut nix, &mut laststart, &mut ixstart,
                     domns, domseglist, alloc, &mut start, nseg);
        }

        // update_ex: compute scores for new segment pairs
        update_ex(
            &mi_w, &trans_w, nseg,
            upper, lower, nseg, segmentrange,
            idom, domns, domseglist, alloc,
            &mut ex, dist, nres1, dist2, nres2,
            dist2sum, dist2sum_stride,
            &mut nix, &mut laststart, &mut ixstart,
            BL, &mut start, &mut lexdone, lseqtl,
            &mut s_beg, s_beg_stride, &mut s_end, s_end_stride,
            dist1sum, dist1sum_stride, minseglen, scoretable,
        );

        // Get initial estimate
        let (ess_init, est_init) = get_ess(ns, &seglist, &mi_w, &ci0_w, nseg,
                                            &ex, &start, nseg,
                                            &mi_w, &trans_w, nseg, lseqtl);

        if ns == 1 {
            // output1: enumerate all singlet candidates above cutoff
            ali_laststart = output1(
                seglist[0], mi_w[seglist[0]] as usize, &ex, cut_w[idom],
                &mut ali_nali, &mut ali_start, ali_laststart, &mut ali_save,
                &trans_w, nseg, idom, &start);

        } else if ns == 2 {
            // output2: enumerate all doublet combinations above cutoff
            ali_laststart = output2(
                seglist[0], seglist[1], &ex, cut_w[idom],
                &mi_w, &ci0_w, nseg, idom, &mut ali_nali, &mut ali_start,
                ali_laststart, &mut ali_save, &trans_w, &start, lseqtl);

        } else {
            // Full branch-and-bound
            let mut pstack = PriorityStack::new();

            // loadstack: initialize priority queue from child alignments
            loadstack(&mut pstack, &ali_start, &ali_nali, &ali_save,
                      &mi_w, &ci0_w, &trans_w, nseg, idom, domns, domseglist, alloc,
                      &ex, &start, lfix, idom1, idom2, ns, &seglist,
                      cut_w[idom], lseqtl, lnul);

            // Memory for degenerate detection
            let mut alimem: HashMap<Vec<i32>, bool> = HashMap::new();

            let mut flag = 5;
            while flag == 5 {
                // dostack: run branch-and-bound search (lprint=false for 10-block)
                let (f, est, ali, ali1) = dostack(
                    &mut pstack, &mi_w, ns, &seglist, cut_w[idom], &ex,
                    &trans_w, nseg, &start, idom, lseqtl,
                    &s_beg, s_beg_stride, &s_end, s_end_stride,
                    segmentrange0, string, outputunit_lines, false);
                flag = f;

                if flag != 5 { break; }

                let est10 = est;
                let mut ali1 = ali1;

                // declump: remove found alignment from stack
                declump(&mut pstack, cut_w[idom], ns, &seglist, &mi_w, &ali,
                        &ex, &start, est10, &trans_w, lseqtl);

                // Check degenerate
                let ali_key: Vec<i32> = seglist.iter()
                    .filter(|&&seg| ali[seg] != NUL)
                    .map(|&seg| ali[seg])
                    .collect();
                if alimem.contains_key(&ali_key) { continue; }
                alimem.insert(ali_key, true);

                // findhomolog: find matching segments in child alignments
                let mut ali2 = vec![NUL; nseg];
                findhomolog(idom1 as i32, &ali1, &ali_nali1, &ali_start1,
                            &ali_save1, &mut ali2, domns, domseglist, alloc, lfix1, nseg);
                findhomolog(idom2 as i32, &ali1, &ali_nali1, &ali_start1,
                            &ali_save1, &mut ali2, domns, domseglist, alloc, lfix1, nseg);

                // perresiduescore: refine alignment at single-residue resolution
                let (est_pr, ali1_pr) = perresiduescore(
                    ns, &seglist, cut_w[idom], lnul, &ali1,
                    upper, lower, segmentrange, dist, nres1, dist2,
                    nres2, dist2sum, ss, nseg, idom, domns, domseglist, alloc,
                    lseqtl, string, &ali2, dist1sum, minseglen,
                    segmentrange0, outputunit_lines, scoretable);
                ali1 = ali1_pr;

                if est_pr <= cut_w[idom] { continue; }

                // saveit: store alignment for parent nodes
                let mut alix = ali.clone();
                for is_ in 0..ns {
                    let iseg = seglist[is_];
                    if ali1[iseg] == NUL && trans_w[ali[iseg] as usize * nseg + iseg] != NUL {
                        // Find NUL candidate
                        for ir in (0..mi_w[iseg] as usize).rev() {
                            if trans_w[ci0_w[ir * nseg + iseg] as usize * nseg + iseg] == NUL {
                                alix[iseg] = mi_w[iseg] - 1;
                                break;
                            }
                        }
                    }
                }

                ali_laststart = saveit(
                    idom, &mut ali_nali, &mut ali_start, ali_laststart,
                    &mut ali_save, &trans_w, nseg, ns, &seglist, &alix);

                // Store perresiduescore result
                if ali_nali1[idom] == 0 {
                    ali_start1[idom] = ali_laststart1;
                }
                ali_nali1[idom] += 1;
                for is_ in 0..ns {
                    let iseg = seglist[is_];
                    ali_laststart1 += 1;
                    if (ali_laststart1 as usize) < ali_save1.len() {
                        ali_save1[ali_laststart1 as usize] = ali1[iseg];
                    }
                }

                if lfirstonly {
                    if idom == 1 { return; }
                    icyc += 1;
                    if icyc >= 5 { break; }
                }
            }
        }
    }
}

/// Run branch-and-bound: pop best, split until unique or exception.
/// Returns (flag, est, ali, ali1) where ali1 has real residue numbers.
fn dostack(
    pstack: &mut PriorityStack, mi: &[i32], ns: usize, seglist: &[usize],
    scorecutoff: i32, ex: &[i32], trans: &[i32], trans_stride: usize,
    start: &[i32], idom: usize, lseqtl: bool,
    s_beg: &[i16], s_beg_stride: usize, s_end: &[i16], s_end_stride: usize,
    segmentrange: &[i32], string: &str,
    outputunit_lines: &mut Vec<String>, lprint: bool,
) -> (i32, i32, Vec<i32>, Vec<i32>) {
    let nseg = mi.len();
    let mut ali = vec![0i32; nseg];
    let mut ali1 = vec![NUL; nseg];

    loop {
        let (flag, est, a) = getnextbest(
            pstack, scorecutoff, ns, seglist, mi, ex, start, trans, lseqtl);
        ali = a;

        match flag {
            1 => return (flag, est, ali, ali1),  // empty stack
            2 => return (flag, est, ali, ali1),  // below cutoff
            3 => {
                // overflow: clearstack, remove bottom 90%
                if pstack.size() <= 1 {
                    return (4, est, ali, ali1);
                }
                let idx_10pct = pstack.size() / 10;
                if idx_10pct < pstack.size() {
                    let k = pstack.states[pstack.size() - idx_10pct - 1].est;
                    pstack.clearstack(k);
                }
                continue;
            }
            4 => continue,  // no candidates
            5 => {
                // unique alignment — convert to real residue numbers
                for is_ in 0..ns {
                    let iseg = seglist[is_];
                    ali1[iseg] = trans[ali[iseg] as usize * trans_stride + iseg];
                }

                if lprint {
                    write_refine(ns, seglist, &ali1, s_beg, s_beg_stride,
                                 s_end, s_end_stride, segmentrange, idom, est,
                                 string, outputunit_lines);
                }

                return (flag, est, ali, ali1);
            }
            _ => continue,
        }
    }
}

/// Write refine output line.
fn write_refine(
    ns: usize, seglist: &[usize], ali1: &[i32],
    s_beg: &[i16], s_beg_stride: usize, s_end: &[i16], s_end_stride: usize,
    segmentrange: &[i32], idom: usize, est: i32,
    string: &str, outputunit_lines: &mut Vec<String>,
) {
    let mut a1 = vec![0i32; ns];
    let mut a2 = vec![0i32; ns];
    let mut b1 = vec![0i32; ns];
    let mut b2 = vec![0i32; ns];

    let s_beg_nrows = s_beg.len() / s_beg_stride;

    for is_ in 0..ns {
        let iseg = seglist[is_];
        let ires = ali1[iseg];
        if ires != NUL {
            let ibeg_key = (ires + 29) as usize;
            let (ibeg, iend) = if ibeg_key < s_beg_nrows {
                (s_beg[ibeg_key * s_beg_stride + iseg] as i32,
                 s_end[ibeg_key * s_end_stride + iseg] as i32)
            } else {
                (0, 0)
            };
            a1[is_] = segmentrange[iseg * 2] + ibeg;
            a2[is_] = segmentrange[iseg * 2 + 1] - iend;
            b1[is_] = ires + ibeg;
            b2[is_] = b1[is_] + a2[is_] - a1[is_];
        } else {
            a1[is_] = NUL;
            a2[is_] = NUL;
            b1[is_] = NUL;
            b2[is_] = NUL;
        }
    }

    let mut parts = String::new();
    parts.push_str(&format!(" refine{}", string));
    parts.push_str(&format!(" {:>10}", idom));
    parts.push_str(&format!(" {:>10}", est));
    parts.push_str(&format!(" {:>10}", ns));
    for is_ in 0..ns {
        parts.push_str(&format!(" {:>10}", a1[is_]));
        parts.push_str(&format!(" {:>10}", a2[is_]));
    }
    for is_ in 0..ns {
        parts.push_str(&format!(" {:>10}", b1[is_]));
        parts.push_str(&format!(" {:>10}", b2[is_]));
    }

    outputunit_lines.push(parts);
}

/// Enumerate all singlet candidates above cutoff.
fn output1(
    iseg: usize, nir: usize, ex: &[i32], cut: i32,
    ali_nali: &mut [i32], ali_start: &mut [i32], mut ali_laststart: i32,
    ali_save: &mut [i32], trans: &[i32], trans_stride: usize,
    idom: usize, start: &[i32],
) -> i32 {
    let seglist = vec![iseg];
    let nseg = trans_stride;
    let mut ali = vec![0i32; nseg.max(iseg + 1)];

    for ir in 0..nir {
        let iwhere = start[iseg * nseg + iseg];
        if iwhere < 0 { continue; }
        let idx = iwhere as usize + ir;
        if idx < ex.len() && ex[idx] > cut {
            ali[iseg] = ir as i32;
            ali_laststart = saveit(
                idom, ali_nali, ali_start, ali_laststart,
                ali_save, trans, trans_stride, 1, &seglist, &ali);
        }
    }
    ali_laststart
}

/// Enumerate all doublet combinations above cutoff.
fn output2(
    aseg: usize, bseg: usize, ex: &[i32], cut: i32,
    mi: &[i32], ci0: &[i32], ci_stride: usize,
    idom: usize, ali_nali: &mut [i32], ali_start: &mut [i32],
    mut ali_laststart: i32, ali_save: &mut [i32],
    trans: &[i32], start: &[i32], lseqtl: bool,
) -> i32 {
    let iseg = aseg.min(bseg);
    let jseg = aseg.max(bseg);
    let nseg = ci_stride;
    let iistart = start[iseg * nseg + iseg];
    let jjstart = start[jseg * nseg + jseg];
    let ijstart = start[iseg * nseg + jseg];
    let seglist = vec![aseg, bseg];
    let jw = mi[jseg];
    let ali_size = aseg.max(bseg) + 1;
    let mut ali = vec![0i32; ali_size];

    for ir in 0..mi[iseg] as usize {
        let ires = ci0[ir * nseg + iseg];
        let iwhere = ijstart + ires * jw;
        for jr in 0..mi[jseg] as usize {
            let jres = ci0[jr * nseg + jseg];
            if lseqtl && trans[jres as usize * nseg + jseg] < trans[ires as usize * nseg + iseg] {
                continue;
            }
            let ii_idx = (iistart + ires) as usize;
            let jj_idx = (jjstart + jres) as usize;
            let ij_idx = (iwhere + jres) as usize;
            if ii_idx >= EXDIM || jj_idx >= EXDIM || ij_idx >= EXDIM {
                continue;
            }
            let mut x = ex[ii_idx] + ex[jj_idx];
            if iwhere >= 0 {
                x += ex[ij_idx];
            }
            if x > cut {
                ali[iseg] = ires;
                ali[jseg] = jres;
                ali_laststart = saveit(
                    idom, ali_nali, ali_start, ali_laststart,
                    ali_save, trans, nseg, 2, &seglist, &ali);
            }
        }
    }
    ali_laststart
}

/// Store alignment for parent nodes. Returns updated ali_laststart.
fn saveit(
    idom: usize, ali_nali: &mut [i32], ali_start: &mut [i32],
    mut ali_laststart: i32, ali_save: &mut [i32],
    trans: &[i32], trans_stride: usize,
    ns: usize, seglist: &[usize], ali: &[i32],
) -> i32 {
    // Check degeneracy
    let ix0 = ali_start[idom] as usize;
    for i in 0..ali_nali[idom] as usize {
        let mut match_found = true;
        for is_ in 0..ns {
            let iseg = seglist[is_];
            let save_idx = ix0 + i * ns + is_ + 1;
            if save_idx < ali_save.len() {
                if trans[ali[iseg] as usize * trans_stride + iseg] != ali_save[save_idx] {
                    match_found = false;
                    break;
                }
            }
        }
        if match_found { return ali_laststart; }
    }

    if ali_nali[idom] == 0 {
        ali_start[idom] = ali_laststart;
    }
    ali_nali[idom] += 1;
    for is_ in 0..ns {
        let iseg = seglist[is_];
        ali_laststart += 1;
        if (ali_laststart as usize) < ali_save.len() {
            ali_save[ali_laststart as usize] = trans[ali[iseg] as usize * trans_stride + iseg];
        }
    }
    ali_laststart
}

/// Build work search space from child node results.
fn shipup(
    ni: &mut Vec<i32>, trans1: &mut Vec<i32>,
    idom: usize, ldom: &[bool], domns: &[i32], domseglist: &[i32], alloc: usize,
    mi: &[i32], ci0: &[i32], trans: &[i32],
    idom1: usize, idom2: usize,
    ali_start: &[i32], ali_nali: &[i32], ali_save: &[i32],
    bl: usize, nres2: usize, lnul: bool,
) {
    let nseg = mi.len();
    *ni = vec![0i32; nseg];
    *trans1 = vec![0i32; trans.len()];

    if idom1 == 0 {
        // No children — load all
        loadall(idom, ni, trans1, mi, ci0, trans, domns, domseglist, alloc, nseg);
    } else {
        if idom1 < ldom.len() && ldom[idom1] {
            loadali(idom1, domns, domseglist, alloc, 0, 0, nres2,
                    ni, trans1, nseg, bl, ali_start, ali_nali, ali_save);
        } else {
            loadall(idom1, ni, trans1, mi, ci0, trans, domns, domseglist, alloc, nseg);
        }
        if idom2 < ldom.len() && ldom[idom2] {
            loadali(idom2, domns, domseglist, alloc, 0, 0, nres2,
                    ni, trans1, nseg, bl, ali_start, ali_nali, ali_save);
        } else {
            loadall(idom2, ni, trans1, mi, ci0, trans, domns, domseglist, alloc, nseg);
        }
    }

    // Add NULs if missing
    if lnul {
        let ns = domns[idom] as usize;
        for is_ in 0..ns {
            let iseg = (domseglist[is_ * alloc + idom] - 1) as usize;
            let mut has_nul = false;
            for ir in 0..ni[iseg] as usize {
                if trans1[ir * nseg + iseg] == NUL {
                    has_nul = true;
                    break;
                }
            }
            if !has_nul {
                trans1[ni[iseg] as usize * nseg + iseg] = NUL;
                ni[iseg] += 1;
            }
        }
    }
}

fn loadall(
    idom: usize, ni: &mut [i32], trans1: &mut [i32],
    mi: &[i32], ci0: &[i32], trans: &[i32],
    domns: &[i32], domseglist: &[i32], alloc: usize, nseg: usize,
) {
    let ns = domns[idom] as usize;
    for is_ in 0..ns {
        let iseg = (domseglist[is_ * alloc + idom] - 1) as usize;
        ni[iseg] = mi[iseg];
        for ir in 0..ni[iseg] as usize {
            trans1[ir * nseg + iseg] = trans[ir * nseg + iseg];
        }
    }
}

fn loadali(
    idom: usize, domns: &[i32], domseglist: &[i32], alloc: usize,
    fuzz_n: i32, fuzz_c: i32, nres2: usize,
    ni: &mut [i32], trans: &mut [i32], nseg: usize,
    bl: usize, ali_start: &[i32], ali_nali: &[i32], ali_save: &[i32],
) {
    let stres: i32 = -29;
    let ns = domns[idom] as usize;

    // Build candidate mask
    let mut lcand: HashMap<usize, std::collections::HashSet<i32>> = HashMap::new();
    for is_ in 0..ns {
        let iseg = (domseglist[is_ * alloc + idom] - 1) as usize;
        lcand.insert(iseg, std::collections::HashSet::new());
    }

    let mut ix = ali_start[idom] as usize;
    for _iali in 0..ali_nali[idom] {
        for is_ in 0..ns {
            ix += 1;
            if ix >= ali_save.len() { break; }
            let ires = ali_save[ix];
            let iseg = (domseglist[is_ * alloc + idom] - 1) as usize;
            if ires != NUL {
                if let Some(set) = lcand.get_mut(&iseg) {
                    let lo = stres.max(ires - fuzz_n);
                    let hi = (nres2 as i32).min(ires + fuzz_c);
                    for i in lo..=hi {
                        set.insert(i);
                    }
                }
            }
        }
    }

    for is_ in 0..ns {
        let iseg = (domseglist[is_ * alloc + idom] - 1) as usize;
        ni[iseg] = 0;
        let set = lcand.get(&iseg);
        let mut ires = stres;
        while ires <= nres2 as i32 {
            if let Some(s) = set {
                if s.contains(&ires) {
                    trans[ni[iseg] as usize * nseg + iseg] = ires;
                    ni[iseg] += 1;
                }
            }
            ires += bl as i32;
        }
    }
}

/// Remap ex[] indices after merging child nodes.
fn compress(
    idom: usize, ni: &[i32], trans1: &[i32],
    mi: &mut [i32], ci0: &mut [i32], trans: &mut [i32],
    ex: &mut [i32], nix: &mut usize, laststart: &mut i32,
    ixstart: &mut [i32], domns: &[i32], domseglist: &[i32],
    alloc: usize, start: &mut [i32], nseg: usize,
) {
    // Build mapping: old ir -> new ir
    let mut irold = vec![0i32; MAXRES0 * nseg];
    let ns = domns[idom] as usize;
    for is_ in 0..ns {
        let iseg = (domseglist[is_ * alloc + idom] - 1) as usize;
        for ir in 0..ni[iseg] as usize {
            let mut ires = 0usize;
            while ires < mi[iseg] as usize && trans[ires * nseg + iseg] != trans1[ir * nseg + iseg] {
                ires += 1;
            }
            if ires >= mi[iseg] as usize { return; }
            irold[ir * nseg + iseg] = ires as i32;
        }
    }

    // Remap singlets and doublets
    for is_ in 0..ns {
        let iseg = (domseglist[is_ * alloc + idom] - 1) as usize;
        let iwhere = start[iseg * nseg + iseg];
        if iwhere >= 0 {
            for ir in 0..ni[iseg].min((EXDIM as i32 - iwhere).max(0)) as usize {
                let old_idx = (iwhere + irold[ir * nseg + iseg]) as usize;
                let new_idx = (iwhere + ir as i32) as usize;
                if old_idx < EXDIM && new_idx < EXDIM {
                    ex[new_idx] = ex[old_idx];
                }
            }
        }

        for js in 0..is_ {
            let jseg = (domseglist[js * alloc + idom] - 1) as usize;
            let aseg = iseg.min(jseg);
            let bseg = iseg.max(jseg);
            if start[aseg * nseg + bseg] >= 0 {
                for ir in 0..ni[aseg] as usize {
                    let iwhere_new = start[aseg * nseg + bseg] + ir as i32 * ni[bseg];
                    let iwhere_old = start[aseg * nseg + bseg] + irold[ir * nseg + aseg] * mi[bseg];
                    for jr in 0..ni[bseg] as usize {
                        let old_idx = (iwhere_old + irold[jr * nseg + bseg]) as usize;
                        let new_idx = (iwhere_new + jr as i32) as usize;
                        if old_idx < EXDIM && new_idx < EXDIM {
                            ex[new_idx] = ex[old_idx];
                        }
                    }
                }
            }
        }
    }

    // Update mi, ci0, trans
    for is_ in 0..ns {
        let iseg = (domseglist[is_ * alloc + idom] - 1) as usize;
        for ir in 0..ni[iseg] as usize {
            trans[ir * nseg + iseg] = trans1[ir * nseg + iseg];
            ci0[ir * nseg + iseg] = ir as i32;
        }
        mi[iseg] = ni[iseg];
    }

    // Left-shift ex chunks
    let mut new_laststart: i32 = 0;
    for ix in 0..*nix {
        let iseg = ixstart[0 * (nseg * nseg) + ix] as usize;
        let jseg = ixstart[1 * (nseg * nseg) + ix] as usize;
        let old_start = start[iseg * nseg + jseg];
        start[iseg * nseg + jseg] = new_laststart;
        start[jseg * nseg + iseg] = new_laststart;
        let count = if iseg == jseg {
            mi[iseg]
        } else {
            mi[iseg] * mi[jseg]
        };
        for i in 0..count {
            let src = (old_start + i) as usize;
            let dst = (new_laststart + i) as usize;
            if src < EXDIM && dst < EXDIM {
                ex[dst] = ex[src];
            }
        }
        new_laststart += count;
    }

    *laststart = new_laststart;
}

/// Initialize priority queue from child alignments.
fn loadstack(
    pstack: &mut PriorityStack,
    ali_start: &[i32], ali_nali: &[i32], ali_save: &[i32],
    mi: &[i32], ci0: &[i32], trans: &[i32], trans_stride: usize,
    idom: usize, domns: &[i32], domseglist: &[i32], alloc: usize,
    ex: &[i32], start: &[i32],
    lfix: &[bool], idom1: usize, idom2: usize,
    ns: usize, seglist: &[usize], cutoff: i32,
    lseqtl: bool, lnul: bool,
) {
    let nseg = mi.len();
    let lfix_stride = alloc;  // lfix indexed as [iseg * alloc + idom]

    // Get fixed segments from children
    let mut slist = Vec::new();
    let mut tlist = Vec::new();
    if idom1 > 0 {
        for is_ in 0..domns[idom1] as usize {
            let iseg_1 = domseglist[is_ * alloc + idom1] as usize;
            if iseg_1 >= 1 {
                let idx = (iseg_1 - 1) * lfix_stride + idom1;
                if idx < lfix.len() && lfix[idx] {
                    slist.push(is_);
                }
            }
        }
    }
    if idom2 > 0 {
        for is_ in 0..domns[idom2] as usize {
            let iseg_1 = domseglist[is_ * alloc + idom2] as usize;
            if iseg_1 >= 1 {
                let idx = (iseg_1 - 1) * lfix_stride + idom2;
                if idx < lfix.len() && lfix[idx] {
                    tlist.push(is_);
                }
            }
        }
    }

    let nali1 = slist.len();
    let nali2 = tlist.len();
    let ia1 = if nali1 > 1 { ali_nali[idom1].max(1) as usize } else { 1 };
    let ia2 = if nali2 > 1 { ali_nali[idom2].max(1) as usize } else { 1 };

    // Build initial search space
    let mut ni_work = vec![0i32; nseg];
    let mut ci_work = vec![0i32; MAXRES0 * nseg];
    for is_ in 0..ns {
        let iseg = seglist[is_];
        ni_work[iseg] = mi[iseg];
        for ir in 0..mi[iseg] as usize {
            ci_work[ir * nseg + iseg] = ci0[ir * nseg + iseg];
        }
    }

    for iali in 0..ia1 {
        // Set fixed segments from child 1
        if nali1 > 1 {
            for k in 0..nali1 {
                let iseg_1 = domseglist[slist[k] * alloc + idom1] as usize;
                let iseg = iseg_1 - 1;
                ni_work[iseg] = 0;
                let save_idx = ali_start[idom1] as usize
                    + iali * domns[idom1] as usize + slist[k] + 1;
                if save_idx < ali_save.len() {
                    let l = ali_save[save_idx];
                    if l != NUL {
                        selectcand(l, iseg, trans, &mut ni_work, &mut ci_work, mi, nseg);
                    }
                    if lnul {
                        selectcand(NUL, iseg, trans, &mut ni_work, &mut ci_work, mi, nseg);
                    }
                }
            }
        }

        for jali in 0..ia2 {
            // Set fixed segments from child 2
            if nali2 > 1 {
                for k in 0..nali2 {
                    let iseg_1 = domseglist[tlist[k] * alloc + idom2] as usize;
                    let iseg = iseg_1 - 1;
                    ni_work[iseg] = 0;
                    let save_idx = ali_start[idom2] as usize
                        + jali * domns[idom2] as usize + tlist[k] + 1;
                    if save_idx < ali_save.len() {
                        let l = ali_save[save_idx];
                        if l != NUL {
                            selectcand(l, iseg, trans, &mut ni_work, &mut ci_work, mi, nseg);
                        }
                        if lnul {
                            selectcand(NUL, iseg, trans, &mut ni_work, &mut ci_work, mi, nseg);
                        }
                    }
                }
            }

            // Push state
            let (_, est) = get_ess(ns, seglist, &ni_work, &ci_work, nseg,
                                    ex, start, nseg, mi, trans, nseg, lseqtl);
            if est > cutoff {
                let mut candidates = HashMap::new();
                for is_ in 0..ns {
                    let seg = seglist[is_];
                    let mut cands = Vec::new();
                    for k in 0..ni_work[seg] as usize {
                        cands.push(ci_work[k * nseg + seg]);
                    }
                    candidates.insert(seg, cands);
                }
                pstack.push(SearchState { est, candidates });
            }
        }
    }
}

fn selectcand(
    l: i32, iseg: usize, trans: &[i32],
    ni: &mut [i32], ci: &mut [i32], mi: &[i32], nseg: usize,
) {
    for ir in 0..mi[iseg] as usize {
        if trans[ir * nseg + iseg] == l {
            ni[iseg] += 1;
            ci[(ni[iseg] - 1) as usize * nseg + iseg] = ir as i32;
            return;
        }
    }
}

/// Find matching segments in child alignments.
fn findhomolog(
    idom1: i32, ali1: &[i32], ali_nali1: &[i32], ali_start1: &[i32],
    ali_save1: &[i32], ali2: &mut [i32],
    domns: &[i32], domseglist: &[i32], alloc: usize,
    lfix1: &[bool], nseg: usize,
) {
    if idom1 <= 0 { return; }
    let idom1 = idom1 as usize;
    let lfix_stride = alloc;

    let ns = domns[idom1] as usize;
    for is_ in 0..ns {
        let iseg = (domseglist[is_ * alloc + idom1] - 1) as usize;
        if iseg < ali2.len() {
            ali2[iseg] = NUL;
        }
    }

    for iali in 0..ali_nali1[idom1] as usize {
        let ix = ali_start1[idom1] as usize + iali * ns;
        let mut found = true;
        for is_ in 0..ns {
            let iseg = (domseglist[is_ * alloc + idom1] - 1) as usize;
            let fix_idx = iseg * lfix_stride + idom1;
            if fix_idx < lfix1.len() && lfix1[fix_idx] {
                let save_idx = ix + is_ + 1;
                if save_idx < ali_save1.len() {
                    let x = ali_save1[save_idx];
                    if x != NUL {
                        if iseg < ali1.len() {
                            if x < ali1[iseg] || x > ali1[iseg] + 9 {
                                found = false;
                                break;
                            }
                        }
                    }
                    if iseg < ali2.len() {
                        ali2[iseg] = x;
                    }
                }
            }
        }
        if found { return; }
    }
}

/// Refine alignment at single-residue resolution.
fn perresiduescore(
    ns: usize, seglist: &[usize], cutoff: i32, lnul: bool, ali1_in: &[i32],
    upper: &[i32], lower: &[i32], segmentrange: &[i32],
    dist: &[i16], nres1: usize, dist2: &[i16], nres2: usize,
    dist2sum: &[i32], ss: &[i32], nseg: usize,
    idom: usize, domns: &[i32], domseglist: &[i32], alloc: usize,
    lseqtl: bool, string: &str, ali2: &[i32],
    dist1sum: &[i32], minseglen: &[i32],
    segmentrange0: &[i32],
    outputunit_lines: &mut Vec<String>, scoretable: &[i32],
) -> (i32, Vec<i32>) {
    let dist2sum_stride = nres2 + 1;
    let dist1sum_stride = nres1 + 1;

    // Initialize fresh scoring state
    let (mut lexdone, mut start) = initlexdonestart(nseg, ss);
    let mut ex = vec![0i32; EXDIM];
    let mut laststart: i32 = 0;
    let mut nix: usize = 0;
    let mut ixstart = vec![0i32; 2 * nseg * nseg];

    let s_beg_stride = nseg;
    let s_end_stride = nseg;
    let s_beg_rows = nres2 + 30;
    let mut s_beg = vec![0i16; s_beg_rows * s_beg_stride];
    let mut s_end = vec![0i16; s_beg_rows * s_end_stride];

    let mut trans_l = vec![0i32; MAXRES0 * nseg];
    let mut ni = vec![0i32; nseg];
    let mut ci = vec![0i32; MAXRES0 * nseg];
    let mut mi_local = vec![0i32; nseg];

    // Build per-residue search space: ±9 residues around 10-block position
    for is_ in 0..ns {
        let iseg = seglist[is_];
        let mut i = 0usize;
        if ali1_in[iseg] != NUL {
            let base = ali1_in[iseg];
            let end = (base + 10).min(nres2 as i32 + 1);
            let mut ires = base;
            while ires < end {
                ci[i * nseg + iseg] = i as i32;
                trans_l[i * nseg + iseg] = ires;
                i += 1;
                ires += 1;
            }
        }
        if lnul {
            ci[i * nseg + iseg] = i as i32;
            trans_l[i * nseg + iseg] = NUL;
            i += 1;
        }
        ni[iseg] = i as i32;
        mi_local[iseg] = i as i32;
    }

    // Overwrite fixed segments with homolog alignment
    for is_ in 0..ns {
        let iseg = seglist[is_];
        if ali2[iseg] != NUL {
            let mut i = 0usize;
            ci[i * nseg + iseg] = i as i32;
            trans_l[i * nseg + iseg] = ali2[iseg];
            i += 1;
            if lnul {
                ci[i * nseg + iseg] = i as i32;
                trans_l[i * nseg + iseg] = NUL;
                i += 1;
            }
            ni[iseg] = i as i32;
            mi_local[iseg] = i as i32;
        }
    }

    // Compute scores with bl=1 (single residue resolution)
    update_ex(
        &mi_local, &trans_l, nseg,
        upper, lower, nseg, segmentrange,
        idom, domns, domseglist, alloc,
        &mut ex, dist, nres1, dist2, nres2,
        dist2sum, dist2sum_stride,
        &mut nix, &mut laststart, &mut ixstart,
        1, &mut start, &mut lexdone, lseqtl,
        &mut s_beg, s_beg_stride, &mut s_end, s_end_stride,
        dist1sum, dist1sum_stride, minseglen, scoretable,
    );

    // Build priority stack — separate first segment candidates
    let mut pstack = PriorityStack::new();
    let (_ess, _est_total) = get_ess(ns, seglist, &ni, &ci, nseg,
                                       &ex, &start, nseg,
                                       &mi_local, &trans_l, nseg, lseqtl);

    let iseg_first = seglist[0];
    let saved_ni = ni[iseg_first];
    for i in 0..mi_local[iseg_first] as usize {
        ni[iseg_first] = 1;
        ci[0 * nseg + iseg_first] = i as i32;
        let (_, est_i) = get_ess(ns, seglist, &ni, &ci, nseg,
                                   &ex, &start, nseg,
                                   &mi_local, &trans_l, nseg, lseqtl);
        if est_i > cutoff {
            let mut candidates = HashMap::new();
            for is_ in 0..ns {
                let seg = seglist[is_];
                let mut cands = Vec::new();
                for k in 0..ni[seg] as usize {
                    cands.push(ci[k * nseg + seg]);
                }
                candidates.insert(seg, cands);
            }
            pstack.push(SearchState { est: est_i, candidates });
        }
    }
    // Restore ni for first seg
    ni[iseg_first] = saved_ni;

    // Initialize ali1 to NUL
    let mut ali1_out = vec![NUL; nseg];

    // Run dostack with lprint=true
    let (flag, est, _, ali1_result) = dostack(
        &mut pstack, &mi_local, ns, seglist, cutoff, &ex,
        &trans_l, nseg, &start, idom, lseqtl,
        &s_beg, s_beg_stride, &s_end, s_end_stride,
        segmentrange0, string, outputunit_lines, true);

    if flag == 5 {
        ali1_out = ali1_result;
    }

    (est, ali1_out)
}
