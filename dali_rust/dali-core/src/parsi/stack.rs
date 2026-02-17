//! Priority queue / stack for PARSI branch-and-bound search.
//!
//! Bit-packed SearchState: 1 bit per candidate per segment (~20 bytes/state
//! vs ~32KB with flat arrays). Matches Fortran's compact representation.

use std::collections::BinaryHeap;
use super::{NUL, INFINIT, MAXRES0};
use super::scoring::{get_ess, get_ess_into, get_estimate};

/// A search state: score estimate + bit-packed candidate masks.
/// Bit layout: segments packed in seglist order, mi[seg] bits per segment.
/// Bit set = candidate present, bit clear = candidate absent.
#[derive(Clone)]
pub struct SearchState {
    pub est: i32,
    pub bits: Vec<u32>,
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
/// Uses BinaryHeap (max-heap) for O(log n) push/pop.
/// Stores bit-packing layout for pack/unpack operations.
pub struct PriorityStack {
    pub states: BinaryHeap<SearchState>,
    pub overflow_limit: usize,
    /// bit_offset[seg] = starting bit position for segment seg's candidates.
    /// Only meaningful for segments in the active seglist.
    bit_offset: Vec<usize>,
    /// Number of u32 words needed per state's bit vector.
    words_per_state: usize,
}

impl PriorityStack {
    /// Create a new priority stack with bit-packing layout for the given domain.
    /// `ns` segments from `seglist`, with `mi[seg]` max candidates per segment.
    pub fn new(ns: usize, seglist: &[usize], mi: &[i32], nseg: usize) -> Self {
        let mut bit_offset = vec![0usize; nseg];
        let mut total_bits = 0usize;
        for is_ in 0..ns {
            let seg = seglist[is_];
            bit_offset[seg] = total_bits;
            total_bits += mi[seg] as usize;
        }
        PriorityStack {
            states: BinaryHeap::new(),
            overflow_limit: 10000,
            bit_offset,
            words_per_state: (total_bits + 31) / 32,
        }
    }

    /// Pack flat ni/ci arrays into a bit-packed SearchState.
    #[inline]
    pub fn pack(&self, est: i32, ni: &[i32], ci: &[i32], nseg: usize,
                ns: usize, seglist: &[usize]) -> SearchState {
        let mut bits = vec![0u32; self.words_per_state];
        for is_ in 0..ns {
            let seg = seglist[is_];
            let offset = self.bit_offset[seg];
            for k in 0..ni[seg] as usize {
                let cand = ci[k * nseg + seg] as usize;
                let bp = offset + cand;
                bits[bp / 32] |= 1 << (bp % 32);
            }
        }
        SearchState { est, bits }
    }

    /// Unpack a bit-packed SearchState into flat ni/ci arrays.
    #[inline]
    pub fn unpack(&self, state: &SearchState, ni: &mut [i32], ci: &mut [i32],
                  nseg: usize, ns: usize, seglist: &[usize], mi: &[i32]) {
        for is_ in 0..ns {
            let seg = seglist[is_];
            ni[seg] = 0;
            let offset = self.bit_offset[seg];
            for ir in 0..mi[seg] as usize {
                let bp = offset + ir;
                if bp / 32 < state.bits.len()
                    && state.bits[bp / 32] & (1 << (bp % 32)) != 0
                {
                    ci[ni[seg] as usize * nseg + seg] = ir as i32;
                    ni[seg] += 1;
                }
            }
        }
    }

    pub fn push(&mut self, state: SearchState) {
        self.states.push(state);
    }

    pub fn pop_best(&mut self) -> Option<SearchState> {
        self.states.pop()
    }

    pub fn peek_best(&self) -> Option<&SearchState> {
        self.states.peek()
    }

    pub fn size(&self) -> usize {
        self.states.len()
    }

    pub fn is_empty(&self) -> bool {
        self.states.is_empty()
    }

    /// Remove all states with estimate <= cutoff.
    pub fn clearstack(&mut self, cutoff: i32) {
        let old = std::mem::take(&mut self.states);
        self.states = old.into_vec().into_iter()
            .filter(|s| s.est > cutoff)
            .collect();
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

    // Scratch buffers — allocated once, reused each iteration
    let mut ni = vec![0i32; nseg];
    let mut ci = vec![0i32; MAXRES0 * nseg];
    let mut ess_scratch = vec![0i32; nseg * nseg];

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

        // Unpack bit-packed state into scratch buffers
        pstack.unpack(&state, &mut ni, &mut ci, nseg, ns, seglist, mi);

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

        // Split (ni/ci are mutable scratch buffers, restored after each copyandput)
        split(est, ns, seglist, ex, scorecutoff, &mut ni, &mut ci, mi, pstack,
              start, trans, lseqtl, &mut ess_scratch);
    }
}

/// Divide search space and push subgroups onto priority stack.
pub fn split(
    _est: i32, ns: usize, seglist: &[usize], ex: &[i32], scorecutoff: i32,
    ni: &mut [i32], ci: &mut [i32], mi: &[i32], pstack: &mut PriorityStack,
    start: &[i32], trans: &[i32], lseqtl: bool,
    ess_scratch: &mut [i32],
) {
    let nseg = mi.len();

    // Check which segments are closed (only 1 candidate)
    // Stack-allocate for typical nseg sizes (≤ 200)
    let mut lclosed_buf = [false; 200];
    let lclosed = &mut lclosed_buf[..nseg];
    for is_ in 0..ns {
        let iseg = seglist[is_];
        lclosed[iseg] = ni[iseg] <= 1;
    }

    // Get segment score contributions (uses pre-allocated ess buffer)
    let est_new = get_ess_into(ns, seglist, ni, ci, nseg, ex, start, nseg,
                                mi, trans, nseg, lseqtl, ess_scratch);
    if est_new <= scorecutoff { return; }

    // Find segment with highest score contribution (getemax)
    let mut xseg: i32 = -1;
    let mut yseg: i32 = -1;
    let mut emaxim = -INFINIT;

    for is_ in 0..ns {
        let iseg = seglist[is_];
        if let Some((e, x, y)) = getemax_update(ess_scratch[iseg * nseg + iseg], emaxim, &lclosed,
                                                   iseg, iseg, xseg, yseg) {
            emaxim = e;
            xseg = x;
            yseg = y;
        }
        for js in 0..is_ {
            let jseg = seglist[js];
            if let Some((e, x, y)) = getemax_update(ess_scratch[iseg * nseg + jseg], emaxim, &lclosed,
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

    // Compute emax per candidate of xseg (stack-allocated)
    let mut emax = [-INFINIT; MAXRES0];
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
                copyandput(cands, xseg, ni, ci, est_new, ns, seglist, ess_scratch,
                           ex, start, mi, trans, lseqtl, scorecutoff, pstack);
            }
        }
    } else {
        if !ci1.is_empty() {
            copyandput(&ci1, xseg, ni, ci, est_new, ns, seglist, ess_scratch,
                       ex, start, mi, trans, lseqtl, scorecutoff, pstack);
        }
    }

    if !ci2.is_empty() {
        copyandput(&ci2, xseg, ni, ci, est_new, ns, seglist, ess_scratch,
                   ex, start, mi, trans, lseqtl, scorecutoff, pstack);
    }
}

#[inline(always)]
fn getemax_update(
    e: i32, emaxim: i32, lclosed: &[bool],
    iseg: usize, jseg: usize, _xseg: i32, _yseg: i32,
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

/// Modify ni/ci in-place for xseg, re-estimate, pack and push.
/// Restores original xseg state before returning.
fn copyandput(
    cand_list: &[i32], xseg: usize, ni: &mut [i32], ci: &mut [i32],
    _est: i32, ns: usize, seglist: &[usize], _ess: &[i32],
    ex: &[i32], start: &[i32], mi: &[i32], trans: &[i32],
    lseqtl: bool, scorecutoff: i32, pstack: &mut PriorityStack,
) {
    if cand_list.is_empty() { return; }
    let nseg = mi.len();

    // Save original xseg state (small: typically 5-40 entries)
    let old_ni_xseg = ni[xseg];
    let old_ci_len = old_ni_xseg as usize;
    let mut old_ci_xseg = Vec::with_capacity(old_ci_len);
    for k in 0..old_ci_len {
        old_ci_xseg.push(ci[k * nseg + xseg]);
    }

    // Overwrite xseg candidates in-place
    ni[xseg] = cand_list.len() as i32;
    for (k, &c) in cand_list.iter().enumerate() {
        ci[k * nseg + xseg] = c;
    }

    // Re-estimate with modified ni/ci
    let f = re_estimate(_est, ns, seglist, _ess, xseg, ex, ni, ci, start,
                        mi, trans, lseqtl);

    if f > scorecutoff {
        // Pack modified state to bits and push
        let state = pstack.pack(f, ni, ci, nseg, ns, seglist);
        pstack.push(state);
    }

    // Restore original xseg state
    ni[xseg] = old_ni_xseg;
    for (k, &c) in old_ci_xseg.iter().enumerate() {
        ci[k * nseg + xseg] = c;
    }
}

#[inline(always)]
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
    let old_states = std::mem::take(&mut pstack.states).into_vec();
    let mut new_states: Vec<SearchState> = Vec::with_capacity(old_states.len());

    // Scratch buffers — allocated once, reused for each state
    let mut ni_tmp = vec![0i32; nseg];
    let mut ci_tmp = vec![0i32; MAXRES0 * nseg];

    for state in old_states {
        let est = state.est;

        if est > singletcutoff {
            new_states.push(state);
            continue;
        }

        // Unpack into scratch buffers
        pstack.unpack(&state, &mut ni_tmp, &mut ci_tmp, nseg, ns, seglist, mi);

        // Remove ali candidates
        let mut lchange = false;
        let mut valid = true;

        for is_ in 0..ns {
            let seg = seglist[is_];
            let m = ali[seg];
            let count = ni_tmp[seg] as usize;
            let mut found = None;
            for k in 0..count {
                if ci_tmp[k * nseg + seg] == m {
                    found = Some(k);
                    break;
                }
            }
            if let Some(pos) = found {
                if trans[m as usize * nseg + seg] != NUL {
                    // Swap with last (order doesn't matter for scoring)
                    let last = count - 1;
                    ci_tmp[pos * nseg + seg] = ci_tmp[last * nseg + seg];
                    ni_tmp[seg] -= 1;
                    lchange = true;
                }
            }
            if ni_tmp[seg] == 0 {
                valid = false;
                break;
            }
        }

        if !valid { continue; }

        let est_val;
        if lchange {
            let (_, e) = get_ess(ns, seglist, &ni_tmp, &ci_tmp, nseg, ex, start, nseg,
                                  mi, trans, nseg, lseqtl);
            est_val = e;
        } else {
            est_val = est;
        }

        if est_val > scorecutoff {
            let new_state = pstack.pack(est_val, &ni_tmp, &ci_tmp, nseg, ns, seglist);
            new_states.push(new_state);
        }
    }

    pstack.states = BinaryHeap::from(new_states);
}
