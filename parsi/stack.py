"""Priority queue / stack for PARSI branch-and-bound search.

Python implementation using sorted list instead of Fortran's bit-packed
linked-list approach. Same semantics, simpler data structures.
"""

import bisect
import numpy as np
from .parsizes import NUL, INFINIT, MAXRES0, MAXSEG, EXDIM
from .scoring import get_ess, get_estimate


class SearchState:
    """A search state: score estimate + active candidates per segment."""
    __slots__ = ['est', 'candidates']

    def __init__(self, est, candidates):
        self.est = est
        self.candidates = candidates  # dict: seg_idx -> list of candidate indices

    def __lt__(self, other):
        return self.est < other.est


class PriorityStack:
    """Score-ordered priority queue of search states.

    Maintains states sorted by estimate (ascending), so pop from end = best.
    """

    def __init__(self):
        self.states = []  # sorted by est ascending
        self.overflow_limit = 10000

    def push(self, state):
        bisect.insort(self.states, state)

    def pop_best(self):
        if not self.states:
            return None
        return self.states.pop()

    def peek_best(self):
        if not self.states:
            return None
        return self.states[-1]

    @property
    def size(self):
        return len(self.states)

    @property
    def empty(self):
        return len(self.states) == 0

    def clearstack(self, cutoff):
        """Remove all states with estimate <= cutoff."""
        # Find first index where est > cutoff
        idx = bisect.bisect_right(
            [s.est for s in self.states], cutoff)
        self.states = self.states[idx:]


def getnextbest(pstack, scorecutoff, ns, seglist, mi, ex, start, trans,
                lseqtl):
    """Pop best state, split until unique alignment or exception.

    Returns (flag, est, ali) where:
        flag=1: empty stack
        flag=2: top score below cutoff
        flag=3: overflow
        flag=4: no candidates
        flag=5: unique alignment in ali
    """
    nseg = len(mi)
    ali = np.full(nseg, 0, dtype=np.int32)
    est = -INFINIT

    while True:
        # Check exceptions
        if pstack.empty:
            return 1, -INFINIT, ali
        best = pstack.peek_best()
        if best.est < scorecutoff:
            return 2, best.est, ali
        if pstack.size > pstack.overflow_limit:
            return 3, best.est, ali

        # Pop best state
        state = pstack.pop_best()
        est = state.est

        # Build ni/ci from state
        ni = np.zeros(nseg, dtype=np.int32)
        ci = np.zeros((MAXRES0, nseg), dtype=np.int32)
        for seg in seglist:
            cands = state.candidates.get(seg, [])
            ni[seg] = len(cands)
            for k, c in enumerate(cands):
                ci[k, seg] = c

        # Check nmin/nmax
        nmin = min(ni[seglist[is_]] for is_ in range(ns))
        nmax = max(ni[seglist[is_]] for is_ in range(ns))

        if nmin == 0:
            continue  # flag=4, skip
        if nmin == 1 and nmax == 1:
            # Unique alignment
            for is_ in range(ns):
                ali[seglist[is_]] = ci[0, seglist[is_]]
            return 5, est, ali

        # Split
        split(est, ns, seglist, ex, scorecutoff, ni, ci, mi, pstack,
              start, trans, lseqtl)


def split(est, ns, seglist, ex, scorecutoff, ni, ci, mi, pstack,
          start, trans, lseqtl):
    """Divide search space and push subgroups onto priority stack."""
    nseg = len(mi)

    # Check which segments are closed (only 1 candidate)
    lclosed = np.zeros(nseg, dtype=bool)
    for is_ in range(ns):
        iseg = seglist[is_]
        lclosed[iseg] = (ni[iseg] <= 1)

    # Get segment score contributions
    ess, est_new = get_ess(ns, seglist, ni, ci, ex, start, mi, trans, lseqtl)
    if est_new <= scorecutoff:
        return

    # Find segment with highest score contribution (getemax)
    xseg = -1
    yseg = -1
    emaxim = -INFINIT

    for is_ in range(ns):
        iseg = seglist[is_]
        result = _getemax_update(ess[iseg, iseg], emaxim, lclosed, iseg, iseg,
                                 xseg, yseg)
        if result is not None:
            emaxim, xseg, yseg = result
        for js in range(is_):
            jseg = seglist[js]
            result = _getemax_update(ess[iseg, jseg], emaxim, lclosed,
                                     iseg, jseg, xseg, yseg)
            if result is not None:
                emaxim, xseg, yseg = result

    if xseg < 0 and yseg < 0:
        return

    # Compute emax per candidate of xseg
    emax = np.full(MAXRES0, -INFINIT, dtype=np.int32)
    if start[xseg, yseg] < 0:
        n1 = 0
    else:
        if xseg == yseg:
            iwhere = start[xseg, xseg]
            for xr in range(ni[xseg]):
                xres = ci[xr, xseg]
                idx = iwhere + xres
                if 0 <= idx < len(ex):
                    emax[xres] = ex[idx]
        elif xseg < yseg:
            ijstart = start[xseg, yseg]
            jw = mi[yseg]
            for xr in range(ni[xseg]):
                e = -INFINIT
                xres = ci[xr, xseg]
                iwhere = ijstart + xres * jw
                for yr in range(ni[yseg]):
                    idx = iwhere + ci[yr, yseg]
                    if 0 <= idx < len(ex):
                        e = max(e, ex[idx])
                emax[xres] = e
        else:  # xseg > yseg
            for xr in range(ni[xseg]):
                emax[ci[xr, xseg]] = -INFINIT
            ijstart = start[yseg, xseg]
            jw = mi[xseg]
            for yr in range(ni[yseg]):
                iwhere = ijstart + ci[yr, yseg] * jw
                for xr in range(ni[xseg]):
                    xres = ci[xr, xseg]
                    idx = iwhere + xres
                    if 0 <= idx < len(ex):
                        emax[xres] = max(emax[xres], ex[idx])
        n1 = ni[xseg]

    # Find min/max of emax values
    if n1 == 0:
        mini = 0
        maxi = 0
    else:
        mini = INFINIT
        maxi = -INFINIT
        bestxres = ci[0, xseg]
        for xr in range(n1):
            xres = ci[xr, xseg]
            e = emax[xres]
            if xr > 1:
                if e > emax[bestxres]:
                    bestxres = xres
            maxi = max(maxi, e)
            mini = min(mini, e)

    # Split 90%-10%
    midpoint = mini + (maxi - mini) * 9 // 10
    ci1 = []
    ci2 = []
    if mini == maxi:
        for xr in range(0, ni[xseg], 2):
            ci1.append(ci[xr, xseg])
        for xr in range(1, ni[xseg], 2):
            ci2.append(ci[xr, xseg])
    else:
        for xr in range(ni[xseg]):
            xres = ci[xr, xseg]
            if emax[xres] <= midpoint:
                ci2.append(xres)
            else:
                ci1.append(xres)

    # Push subgroups
    if lseqtl and n1 > 0:
        # Split better half into left/nul/right groups
        ci3 = []  # left of bestxres
        ci4 = []  # right of bestxres (including bestxres)
        ci6 = []  # nul candidates
        for xres in ci1:
            pr = trans[xres, xseg]
            if pr == NUL:
                ci6.append(xres)
            elif xres < bestxres:
                ci3.append(xres)
            else:
                ci4.append(xres)
        for cands in [ci6, ci3, ci4]:
            if cands:
                _copyandput(cands, xseg, ni, ci, est_new, ns, seglist, ess,
                            ex, start, mi, trans, lseqtl, scorecutoff, pstack)
    else:
        if ci1:
            _copyandput(ci1, xseg, ni, ci, est_new, ns, seglist, ess,
                        ex, start, mi, trans, lseqtl, scorecutoff, pstack)

    # Push lower half
    if ci2:
        _copyandput(ci2, xseg, ni, ci, est_new, ns, seglist, ess,
                    ex, start, mi, trans, lseqtl, scorecutoff, pstack)


def _getemax_update(e, emaxim, lclosed, iseg, jseg, xseg, yseg):
    """Update xseg/yseg if e > emaxim and at least one segment is open.

    Uses gfortran RAND quirk: ran(1234567) always returns 0.0878... < 0.5,
    so when both are open, always picks xseg=jseg, yseg=iseg.
    """
    if e > emaxim:
        if lclosed[iseg] and not lclosed[jseg]:
            return e, jseg, iseg
        elif lclosed[jseg] and not lclosed[iseg]:
            return e, iseg, jseg
        elif not lclosed[iseg] and not lclosed[jseg]:
            # ran(seed) = 0.0878... < 0.5, so else branch
            return e, jseg, iseg
    return None


def _copyandput(cand_list, xseg, ni, ci, est, ns, seglist, ess,
                ex, start, mi, trans, lseqtl, scorecutoff, pstack):
    """Set xseg candidates to cand_list, re-estimate, and push if above cutoff."""
    if not cand_list:
        return

    nseg = len(mi)
    # Save original
    old_ni_xseg = ni[xseg]
    old_ci_xseg = ci[:old_ni_xseg, xseg].copy()

    # Set xseg candidates
    ni[xseg] = len(cand_list)
    for k, c in enumerate(cand_list):
        ci[k, xseg] = c

    # Re-estimate
    f = re_estimate(est, ns, seglist, ess, xseg, ex, ni, ci, start, mi,
                    trans, lseqtl)

    if f > scorecutoff:
        # Build state
        candidates = {}
        for is_ in range(ns):
            seg = seglist[is_]
            candidates[seg] = [ci[k, seg] for k in range(ni[seg])]
        pstack.push(SearchState(f, candidates))

    # Restore original
    ni[xseg] = old_ni_xseg
    ci[:old_ni_xseg, xseg] = old_ci_xseg


def re_estimate(est, ns, seglist, ess, xseg, ex, ni, ci, start, mi,
                trans, lseqtl):
    """Re-estimate score after modifying xseg candidates."""
    nseg = len(mi)
    new_ess = np.zeros((nseg, nseg), dtype=np.int32)
    new_est = 0

    for is_ in range(ns):
        iseg = seglist[is_]
        if iseg == xseg or True:  # Recalculate all involving xseg
            pass
        e = get_estimate(iseg, iseg, ni, ci, ex, start, mi, trans, lseqtl)
        new_ess[iseg, iseg] = e
        new_est += e
        for js in range(is_):
            jseg = seglist[js]
            e = get_estimate(iseg, jseg, ni, ci, ex, start, mi, trans, lseqtl)
            e = e + e
            new_ess[iseg, jseg] = e
            new_ess[jseg, iseg] = e
            new_est += e

    return new_est


def declump(pstack, scorecutoff, ns, seglist, mi, ali, ex, start,
            singletcutoff, trans, lseqtl):
    """Remove found alignment (ali) from all stack states.

    For each state, remove the found candidate for each segment.
    Re-estimate scores for modified states.
    """
    nseg = len(mi)
    new_states = []

    for state in pstack.states:
        est = state.est

        # Only declump states below singletcutoff
        if est > singletcutoff:
            new_states.append(state)
            continue

        # Remove ali candidates from state
        lchange = False
        new_candidates = {}
        valid = True
        for is_ in range(ns):
            seg = seglist[is_]
            cands = list(state.candidates.get(seg, []))
            m = ali[seg]
            # Remove candidate m if it maps to a non-NUL residue
            if m in cands and trans[m, seg] != NUL:
                cands.remove(m)
                lchange = True
            if not cands:
                valid = False
                break
            new_candidates[seg] = cands

        if not valid:
            continue

        if lchange:
            # Re-estimate
            ni_tmp = np.zeros(nseg, dtype=np.int32)
            ci_tmp = np.zeros((MAXRES0, nseg), dtype=np.int32)
            for seg, cands in new_candidates.items():
                ni_tmp[seg] = len(cands)
                for k, c in enumerate(cands):
                    ci_tmp[k, seg] = c
            _, est = get_ess(ns, seglist, ni_tmp, ci_tmp, ex, start, mi,
                             trans, lseqtl)

        if est > scorecutoff:
            new_states.append(SearchState(est, new_candidates))

    pstack.states = sorted(new_states)
