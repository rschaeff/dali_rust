//! Priority queue / stack for PARSI branch-and-bound search.
//!
//! Python-style sorted Vec instead of Fortran's bit-packed linked-list.

use std::collections::HashMap;
use super::{NUL, INFINIT, MAXRES0};
use super::scoring::{get_ess, get_estimate};

/// A search state: score estimate + active candidates per segment.
#[derive(Clone)]
pub struct SearchState {
    pub est: i32,
    pub candidates: HashMap<usize, Vec<i32>>,  // seg_idx -> candidate indices
}

impl PartialEq for SearchState {
    fn eq(&self, other: &Self) -> bool { self.est == other.est }
}
impl Eq for SearchState {}
impl PartialOrd for SearchState {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for SearchState {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.est.cmp(&other.est)
    }
}

/// Score-ordered priority queue of search states.
/// Maintains states sorted by estimate (ascending), so pop from end = best.
pub struct PriorityStack {
    pub states: Vec<SearchState>,
    pub overflow_limit: usize,
}

impl PriorityStack {
    pub fn new() -> Self {
        PriorityStack {
            states: Vec::new(),
            overflow_limit: 10000,
        }
    }

    pub fn push(&mut self, state: SearchState) {
        let pos = self.states.partition_point(|s| s.est <= state.est);
        self.states.insert(pos, state);
    }

    pub fn pop_best(&mut self) -> Option<SearchState> {
        if self.states.is_empty() { None } else { Some(self.states.pop().unwrap()) }
    }

    pub fn peek_best(&self) -> Option<&SearchState> {
        self.states.last()
    }

    pub fn size(&self) -> usize {
        self.states.len()
    }

    pub fn is_empty(&self) -> bool {
        self.states.is_empty()
    }

    /// Remove all states with estimate <= cutoff.
    pub fn clearstack(&mut self, cutoff: i32) {
        let idx = self.states.partition_point(|s| s.est <= cutoff);
        self.states = self.states[idx..].to_vec();
    }
}

/// Pop best state, split until unique alignment or exception.
///
/// Returns (flag, est, ali) where:
///   flag=1: empty stack
///   flag=2: top score below cutoff
///   flag=3: overflow
///   flag=4: no candidates
///   flag=5: unique alignment in ali
pub fn getnextbest(
    pstack: &mut PriorityStack, scorecutoff: i32,
    ns: usize, seglist: &[usize], mi: &[i32],
    ex: &[i32], start: &[i32], trans: &[i32], lseqtl: bool,
) -> (i32, i32, Vec<i32>) {
    let nseg = mi.len();
    let mut ali = vec![0i32; nseg];
    let mut est = -INFINIT;

    loop {
        if pstack.is_empty() {
            return (1, -INFINIT, ali);
        }
        let best = pstack.peek_best().unwrap();
        if best.est < scorecutoff {
            return (2, best.est, ali);
        }
        if pstack.size() > pstack.overflow_limit {
            return (3, best.est, ali);
        }

        let state = pstack.pop_best().unwrap();
        est = state.est;

        // Build ni/ci from state
        let mut ni = vec![0i32; nseg];
        let mut ci = vec![0i32; MAXRES0 * nseg];
        for &seg in &seglist[..ns] {
            let cands = state.candidates.get(&seg).map(|v| v.as_slice()).unwrap_or(&[]);
            ni[seg] = cands.len() as i32;
            for (k, &c) in cands.iter().enumerate() {
                ci[k * nseg + seg] = c;
            }
        }

        // Check nmin/nmax
        let mut nmin = i32::MAX;
        let mut nmax = i32::MIN;
        for is_ in 0..ns {
            let seg = seglist[is_];
            nmin = nmin.min(ni[seg]);
            nmax = nmax.max(ni[seg]);
        }

        if nmin == 0 {
            continue;  // skip, no candidates
        }
        if nmin == 1 && nmax == 1 {
            // Unique alignment
            for is_ in 0..ns {
                ali[seglist[is_]] = ci[0 * nseg + seglist[is_]];
            }
            return (5, est, ali);
        }

        // Split
        split(est, ns, seglist, ex, scorecutoff, &ni, &ci, mi, pstack,
              start, trans, lseqtl);
    }
}

/// Divide search space and push subgroups onto priority stack.
pub fn split(
    est: i32, ns: usize, seglist: &[usize], ex: &[i32], scorecutoff: i32,
    ni: &[i32], ci: &[i32], mi: &[i32], pstack: &mut PriorityStack,
    start: &[i32], trans: &[i32], lseqtl: bool,
) {
    let nseg = mi.len();

    // Check which segments are closed (only 1 candidate)
    let mut lclosed = vec![false; nseg];
    for is_ in 0..ns {
        let iseg = seglist[is_];
        lclosed[iseg] = ni[iseg] <= 1;
    }

    // Get segment score contributions
    let (ess, est_new) = get_ess(ns, seglist, ni, ci, nseg, ex, start, nseg,
                                  mi, trans, nseg, lseqtl);
    if est_new <= scorecutoff { return; }

    // Find segment with highest score contribution (getemax)
    let mut xseg: i32 = -1;
    let mut yseg: i32 = -1;
    let mut emaxim = -INFINIT;

    for is_ in 0..ns {
        let iseg = seglist[is_];
        if let Some((e, x, y)) = getemax_update(ess[iseg * nseg + iseg], emaxim, &lclosed,
                                                   iseg, iseg, xseg, yseg) {
            emaxim = e;
            xseg = x;
            yseg = y;
        }
        for js in 0..is_ {
            let jseg = seglist[js];
            if let Some((e, x, y)) = getemax_update(ess[iseg * nseg + jseg], emaxim, &lclosed,
                                                       iseg, jseg, xseg, yseg) {
                emaxim = e;
                xseg = x;
                yseg = y;
            }
        }
    }

    if xseg < 0 && yseg < 0 { return; }
    let xseg = xseg as usize;
    let yseg = yseg as usize;

    // Compute emax per candidate of xseg
    let mut emax = vec![-INFINIT; MAXRES0];
    let n1;
    if start[xseg * nseg + yseg] < 0 {
        n1 = 0;
    } else {
        if xseg == yseg {
            let iwhere = start[xseg * nseg + xseg];
            for xr in 0..ni[xseg] as usize {
                let xres = ci[xr * nseg + xseg] as usize;
                let idx = iwhere as usize + xres;
                if idx < ex.len() {
                    emax[xres] = ex[idx];
                }
            }
        } else if xseg < yseg {
            let ijstart = start[xseg * nseg + yseg];
            let jw = mi[yseg];
            for xr in 0..ni[xseg] as usize {
                let mut e = -INFINIT;
                let xres = ci[xr * nseg + xseg] as usize;
                let iwhere = ijstart + xres as i32 * jw;
                for yr in 0..ni[yseg] as usize {
                    let idx = (iwhere + ci[yr * nseg + yseg]) as usize;
                    if idx < ex.len() {
                        e = e.max(ex[idx]);
                    }
                }
                emax[xres] = e;
            }
        } else {
            // xseg > yseg
            for xr in 0..ni[xseg] as usize {
                emax[ci[xr * nseg + xseg] as usize] = -INFINIT;
            }
            let ijstart = start[yseg * nseg + xseg];
            let jw = mi[xseg];
            for yr in 0..ni[yseg] as usize {
                let iwhere = ijstart + ci[yr * nseg + yseg] as i32 * jw;
                for xr in 0..ni[xseg] as usize {
                    let xres = ci[xr * nseg + xseg] as usize;
                    let idx = (iwhere + xres as i32) as usize;
                    if idx < ex.len() {
                        emax[xres] = emax[xres].max(ex[idx]);
                    }
                }
            }
        }
        n1 = ni[xseg] as usize;
    }

    // Find min/max of emax values
    let (mini, maxi, bestxres);
    if n1 == 0 {
        mini = 0;
        maxi = 0;
        bestxres = 0;
    } else {
        let mut mn = INFINIT;
        let mut mx = -INFINIT;
        let mut bx = ci[0 * nseg + xseg] as usize;
        for xr in 0..n1 {
            let xres = ci[xr * nseg + xseg] as usize;
            let e = emax[xres];
            if xr > 1 {
                if e > emax[bx] {
                    bx = xres;
                }
            }
            mx = mx.max(e);
            mn = mn.min(e);
        }
        mini = mn;
        maxi = mx;
        bestxres = bx;
    }

    // Split 90%-10%
    let midpoint = mini + (maxi - mini) * 9 / 10;
    let mut ci1 = Vec::new();
    let mut ci2 = Vec::new();
    if mini == maxi {
        for xr in (0..ni[xseg] as usize).step_by(2) {
            ci1.push(ci[xr * nseg + xseg]);
        }
        for xr in (1..ni[xseg] as usize).step_by(2) {
            ci2.push(ci[xr * nseg + xseg]);
        }
    } else {
        for xr in 0..ni[xseg] as usize {
            let xres = ci[xr * nseg + xseg] as usize;
            if emax[xres] <= midpoint {
                ci2.push(xres as i32);
            } else {
                ci1.push(xres as i32);
            }
        }
    }

    // Push subgroups
    if lseqtl && n1 > 0 {
        // Split better half into left/nul/right groups
        let mut ci3 = Vec::new();
        let mut ci4 = Vec::new();
        let mut ci6 = Vec::new();
        for &xres in &ci1 {
            let pr = trans[xres as usize * nseg + xseg];
            if pr == NUL {
                ci6.push(xres);
            } else if (xres as usize) < bestxres {
                ci3.push(xres);
            } else {
                ci4.push(xres);
            }
        }
        for cands in [&ci6, &ci3, &ci4] {
            if !cands.is_empty() {
                copyandput(cands, xseg, ni, ci, est_new, ns, seglist, &ess,
                           ex, start, mi, trans, lseqtl, scorecutoff, pstack);
            }
        }
    } else {
        if !ci1.is_empty() {
            copyandput(&ci1, xseg, ni, ci, est_new, ns, seglist, &ess,
                       ex, start, mi, trans, lseqtl, scorecutoff, pstack);
        }
    }

    if !ci2.is_empty() {
        copyandput(&ci2, xseg, ni, ci, est_new, ns, seglist, &ess,
                   ex, start, mi, trans, lseqtl, scorecutoff, pstack);
    }
}

fn getemax_update(
    e: i32, emaxim: i32, lclosed: &[bool],
    iseg: usize, jseg: usize, xseg: i32, yseg: i32,
) -> Option<(i32, i32, i32)> {
    if e > emaxim {
        if lclosed[iseg] && !lclosed[jseg] {
            return Some((e, jseg as i32, iseg as i32));
        } else if lclosed[jseg] && !lclosed[iseg] {
            return Some((e, iseg as i32, jseg as i32));
        } else if !lclosed[iseg] && !lclosed[jseg] {
            // ran(seed) = 0.0878... < 0.5, so always picks jseg, iseg
            return Some((e, jseg as i32, iseg as i32));
        }
    }
    None
}

fn copyandput(
    cand_list: &[i32], xseg: usize, ni: &[i32], ci: &[i32],
    est: i32, ns: usize, seglist: &[usize], ess: &[i32],
    ex: &[i32], start: &[i32], mi: &[i32], trans: &[i32],
    lseqtl: bool, scorecutoff: i32, pstack: &mut PriorityStack,
) {
    if cand_list.is_empty() { return; }
    let nseg = mi.len();

    // Save original
    let old_ni_xseg = ni[xseg];
    let mut old_ci_xseg = Vec::with_capacity(old_ni_xseg as usize);
    for k in 0..old_ni_xseg as usize {
        old_ci_xseg.push(ci[k * nseg + xseg]);
    }

    // Create mutable copies
    let mut ni_w = ni.to_vec();
    let mut ci_w = ci.to_vec();

    // Set xseg candidates
    ni_w[xseg] = cand_list.len() as i32;
    for (k, &c) in cand_list.iter().enumerate() {
        ci_w[k * nseg + xseg] = c;
    }

    // Re-estimate
    let f = re_estimate(est, ns, seglist, ess, xseg, ex, &ni_w, &ci_w, start,
                        mi, trans, lseqtl);

    if f > scorecutoff {
        let mut candidates = HashMap::new();
        for is_ in 0..ns {
            let seg = seglist[is_];
            let mut cands = Vec::new();
            for k in 0..ni_w[seg] as usize {
                cands.push(ci_w[k * nseg + seg]);
            }
            candidates.insert(seg, cands);
        }
        pstack.push(SearchState { est: f, candidates });
    }
}

fn re_estimate(
    _est: i32, ns: usize, seglist: &[usize], _ess: &[i32], _xseg: usize,
    ex: &[i32], ni: &[i32], ci: &[i32], start: &[i32],
    mi: &[i32], trans: &[i32], lseqtl: bool,
) -> i32 {
    let nseg = mi.len();
    let mut new_est: i64 = 0;
    for is_ in 0..ns {
        let iseg = seglist[is_];
        let e = get_estimate(iseg, iseg, ni, ci, nseg, ex, start, nseg,
                              mi, trans, nseg, lseqtl);
        new_est += e as i64;
        for js in 0..is_ {
            let jseg = seglist[js];
            let e = get_estimate(iseg, jseg, ni, ci, nseg, ex, start, nseg,
                                  mi, trans, nseg, lseqtl);
            new_est += (e + e) as i64;
        }
    }
    new_est.min(i32::MAX as i64) as i32
}

/// Remove found alignment (ali) from all stack states.
pub fn declump(
    pstack: &mut PriorityStack, scorecutoff: i32,
    ns: usize, seglist: &[usize], mi: &[i32], ali: &[i32],
    ex: &[i32], start: &[i32], singletcutoff: i32,
    trans: &[i32], lseqtl: bool,
) {
    let nseg = mi.len();
    let mut new_states = Vec::new();

    for state in &pstack.states {
        let est = state.est;

        if est > singletcutoff {
            new_states.push(state.clone());
            continue;
        }

        // Remove ali candidates from state
        let mut lchange = false;
        let mut new_candidates = HashMap::new();
        let mut valid = true;
        for is_ in 0..ns {
            let seg = seglist[is_];
            let mut cands: Vec<i32> = state.candidates.get(&seg)
                .map(|v| v.clone())
                .unwrap_or_default();
            let m = ali[seg];
            // Remove candidate m if it maps to a non-NUL residue
            if let Some(pos) = cands.iter().position(|&c| c == m) {
                if trans[m as usize * nseg + seg] != NUL {
                    cands.remove(pos);
                    lchange = true;
                }
            }
            if cands.is_empty() {
                valid = false;
                break;
            }
            new_candidates.insert(seg, cands);
        }

        if !valid { continue; }

        let est_val;
        if lchange {
            // Re-estimate
            let mut ni_tmp = vec![0i32; nseg];
            let mut ci_tmp = vec![0i32; MAXRES0 * nseg];
            for (&seg, cands) in &new_candidates {
                ni_tmp[seg] = cands.len() as i32;
                for (k, &c) in cands.iter().enumerate() {
                    ci_tmp[k * nseg + seg] = c;
                }
            }
            let (_, e) = get_ess(ns, seglist, &ni_tmp, &ci_tmp, nseg, ex, start, nseg,
                                  mi, trans, nseg, lseqtl);
            est_val = e;
        } else {
            est_val = est;
        }

        if est_val > scorecutoff {
            new_states.push(SearchState { est: est_val, candidates: new_candidates });
        }
    }

    new_states.sort();
    pstack.states = new_states;
}
