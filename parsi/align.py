"""Core alignment module for PARSI — align() orchestrator + helpers."""

import numpy as np
from .parsizes import (NUL, INFINIT, EXDIM, MAXRES0, MAXSEG, MAXDOM,
                       MAXSTACK, MAXSTACK2, BL, BOXDIM)
from .scoring import (get_ess, get_estimate, update_ex, initlexdonestart,
                      singletex, doubletex, segsegscore)
from .stack import PriorityStack, SearchState, getnextbest, declump


def align(prot1, prot2_nres, prot2_secstr, prot2_nseg, prot2_segment,
          dist, nres1, dist2, nres2, dist2sum, dist1sum,
          ss, upper, lower, segmentrange, segmentrange0,
          ndom, node_child, domns, domseglist, ldom, lfix, lfix1,
          cut, minseglen, ngap, checkrange, checkx,
          mi, ci0, trans, scoretable, weight,
          lseqtl, lnul, lfirstonly, string, outputunit_lines):
    """Main PARSI alignment — bottom-up over domain tree.

    Args:
        outputunit_lines: list to append output refine lines to.

    All arrays use 0-based indexing for Python arrays. Segment indices in
    domseglist are 1-based (matching Fortran convention).
    """
    nseg = prot1.nseg

    # Initialize ex array and tracking
    lexdone, start = initlexdonestart(nseg, ss)
    ex = np.zeros(EXDIM, dtype=np.int32)
    nix = 0
    laststart = 0
    ixstart = np.zeros((2, nseg * nseg), dtype=np.int32)

    # Per-residue score trimming arrays
    s_beg = np.zeros((nres2 + 30, nseg), dtype=np.int16)  # indexed by transires+29
    s_end = np.zeros((nres2 + 30, nseg), dtype=np.int16)

    # Ali save arrays for child results
    ali_save = np.zeros(100001, dtype=np.int32)
    ali_start = np.zeros(ndom + 100, dtype=np.int32)
    ali_nali = np.zeros(ndom + 100, dtype=np.int32)
    ali_laststart = 0

    # Ali save for per-residue (findhomolog)
    ali_save1 = np.zeros(100001, dtype=np.int32)
    ali_start1 = np.zeros(ndom + 100, dtype=np.int32)
    ali_nali1 = np.zeros(ndom + 100, dtype=np.int32)
    ali_laststart1 = 0

    # Working copy of trans for shipup
    trans1 = np.zeros_like(trans)
    ni = np.zeros(nseg, dtype=np.int32)

    # Cut working copy
    cut_w = cut.copy()

    # Process domain tree bottom-up
    _debug = False  # TEMP DEBUG
    for idom in range(ndom, 0, -1):
        icyc = 0
        if not ldom[idom]:
            continue
        if _debug:
            print(f'  [align] idom={idom}, ns={domns[idom]}', flush=True)

        idom1 = node_child[0, idom]
        idom2 = node_child[1, idom]
        ns = domns[idom]
        seglist = np.array([domseglist[i, idom] for i in range(ns)], dtype=np.int32)
        # Convert to 0-based for array access
        seglist_0 = seglist - 1

        # shipup: build work search space from child node results
        ni, trans1 = shipup(ni, trans1, idom, ldom, domns, domseglist,
                            mi, ci0, trans, idom1, idom2,
                            ali_start, ali_nali, ali_save, BL, nres2, lnul)

        # compress: remap ex[] indices after merging child nodes
        if idom1 > 0:
            nix, laststart = compress(idom1, ni, trans1, mi, ci0, trans,
                                       ex, nix, laststart, ixstart,
                                       domns, domseglist, start)
        if idom2 > 0:
            nix, laststart = compress(idom2, ni, trans1, mi, ci0, trans,
                                       ex, nix, laststart, ixstart,
                                       domns, domseglist, start)

        # update_ex: compute scores for new segment pairs
        nix, laststart = update_ex(
            mi, trans, upper, lower, segmentrange, idom, domns,
            domseglist, ex, dist, nres1, dist2, nres2, dist2sum,
            nix, laststart, ixstart, BL, start, lexdone, lseqtl,
            s_beg, s_end, dist1sum, minseglen, scoretable)

        # Get initial estimate
        ess_init, est_init = get_ess(ns, seglist_0, mi, ci0, ex, start,
                                     mi, trans, lseqtl)
        firstest = float(max(est_init, 1))
        if _debug:
            print(f'  [align] idom={idom}: est_init={est_init}, cut={cut_w[idom]}, laststart={laststart}', flush=True)
            print(f'    mi={[mi[s] for s in seglist_0]}, ali_nali={[(i, ali_nali[i]) for i in range(1, ndom+1) if ldom[i]]}', flush=True)

        if ns == 1:
            # output1: enumerate all singlet candidates above cutoff
            iseg_0 = seglist_0[0]
            ali_laststart = output1(
                iseg_0, mi[iseg_0], ex, cut_w[idom],
                ali_nali, ali_start, ali_laststart, ali_save,
                trans, idom, start)

        elif ns == 2:
            # output2: enumerate all doublet combinations above cutoff
            ali_laststart = output2(
                seglist_0[0], seglist_0[1], ex, cut_w[idom],
                mi, ci0, idom, ali_nali, ali_start, ali_laststart,
                ali_save, trans, start, lseqtl)

        else:
            # Full branch-and-bound
            pstack = PriorityStack()

            # loadstack: initialize priority queue from child alignments
            loadstack(pstack, ali_start, ali_nali, ali_save,
                      mi, ci0, trans, idom, domns, domseglist, ex,
                      start, lfix, idom1, idom2, ns, seglist_0,
                      cut_w[idom], lseqtl, lnul)

            # Memory for degenerate detection
            alimem = {}
            nmem = 0

            flag = 5
            while flag == 5:
                # dostack: run branch-and-bound search (lprint=False for 10-block)
                flag, est, ali, ali1 = dostack(
                    pstack, mi, ns, seglist_0, cut_w[idom], ex,
                    trans, start, idom, lseqtl, s_beg, s_end,
                    segmentrange0, string, outputunit_lines,
                    lprint=False)
                if _debug:
                    print(f'  [align] idom={idom}: dostack flag={flag}, est={est}', flush=True)

                if flag != 5:
                    break

                est10 = est

                if _debug and idom <= 2:
                    ali1_str = ','.join(str(ali1[s]) for s in seglist_0)
                    print(f'  [align] idom={idom} cycle: 10-blk est={est10}, ali1=[{ali1_str}]', flush=True)
                    print(f'  [align] before declump: pstack_size={pstack.size}', flush=True)

                # declump: remove found alignment from stack
                declump(pstack, cut_w[idom], ns, seglist_0, mi, ali,
                        ex, start, est10, trans, lseqtl)

                if _debug and idom <= 2:
                    print(f'  [align] after declump: pstack_size={pstack.size}', flush=True)
                    # Show top 10 states
                    top_n = pstack.states[-10:] if len(pstack.states) >= 10 else pstack.states[:]
                    for idx, s in enumerate(reversed(top_n)):
                        cstr = ', '.join(f's{k}:{len(v)}' for k, v in sorted(s.candidates.items()))
                        print(f'    [{idx}] est={s.est}  cands: {cstr}', flush=True)

                # Check degenerate
                ali_key = tuple(ali[seg] for seg in seglist_0
                                if ali[seg] != NUL)
                if ali_key in alimem:
                    continue
                alimem[ali_key] = True
                nmem += 1

                # findhomolog: find matching segments in child alignments
                ali2 = np.full(nseg, NUL, dtype=np.int32)
                findhomolog(idom1, ali1, ali_nali1, ali_start1,
                            ali_save1, ali2, domns, domseglist, lfix1)
                findhomolog(idom2, ali1, ali_nali1, ali_start1,
                            ali_save1, ali2, domns, domseglist, lfix1)

                # perresiduescore: refine alignment at single-residue resolution
                est, ali1 = perresiduescore(
                    ns, seglist_0, cut_w[idom], lnul, ali1,
                    upper, lower, segmentrange, dist, nres1, dist2,
                    nres2, dist2sum, ss, nseg, idom, domns, domseglist,
                    lseqtl, string, ali2, dist1sum, minseglen,
                    segmentrange0, outputunit_lines, scoretable)

                if _debug:
                    print(f'  [align] idom={idom}: perres est={est}, cut={cut_w[idom]}', flush=True)
                if est <= cut_w[idom]:
                    continue

                # saveit: store alignment for parent nodes
                # Build alix — handle transversions
                alix = ali.copy()
                for is_ in range(ns):
                    iseg = seglist_0[is_]
                    if ali1[iseg] == NUL and trans[ali[iseg], iseg] != NUL:
                        # Find NUL candidate
                        for ir in range(mi[iseg] - 1, -1, -1):
                            if trans[ci0[ir, iseg], iseg] == NUL:
                                alix[iseg] = mi[iseg] - 1  # 0-based last
                                break

                ali_laststart = saveit(
                    idom, ali_nali, ali_start, ali_laststart,
                    ali_save, trans, ns, seglist_0, alix)

                # Store perresiduescore result
                if ali_nali1[idom] == 0:
                    ali_start1[idom] = ali_laststart1
                ali_nali1[idom] += 1
                for is_ in range(ns):
                    iseg = seglist_0[is_]
                    ali_laststart1 += 1
                    if ali_laststart1 < len(ali_save1):
                        ali_save1[ali_laststart1] = ali1[iseg]

                if lfirstonly:
                    if idom == 1:
                        return
                    icyc += 1
                    if icyc >= 5:
                        break


def dostack(pstack, mi, ns, seglist, scorecutoff, ex, trans, start,
            idom, lseqtl, s_beg, s_end, segmentrange, string,
            outputunit_lines, lprint=True):
    """Run branch-and-bound: pop best, split until unique or exception.

    Returns (flag, est, ali, ali1) where ali1 has real residue numbers.
    lprint: if True, write refine output lines (False for 10-block search).
    """
    nseg = len(mi)
    ali = np.zeros(nseg, dtype=np.int32)
    ali1 = np.full(nseg, NUL, dtype=np.int32)

    while True:
        flag, est, ali = getnextbest(
            pstack, scorecutoff, ns, seglist, mi, ex, start,
            trans, lseqtl)

        if flag == 1:  # empty stack
            return flag, est, ali, ali1
        elif flag == 2:  # top score below cutoff
            return flag, est, ali, ali1
        elif flag == 3:  # overflow
            # clearstack: remove bottom 90%
            if pstack.size <= 1:
                return 4, est, ali, ali1
            # Find score at 10th percentile
            idx_10pct = pstack.size // 10
            if idx_10pct < pstack.size:
                k = pstack.states[-(idx_10pct + 1)].est
                pstack.clearstack(k)
            continue
        elif flag == 4:  # no candidates
            continue
        elif flag == 5:  # unique alignment
            # Convert to real residue numbers
            for is_ in range(ns):
                iseg = seglist[is_]
                ali1[iseg] = trans[ali[iseg], iseg]

            # Write output only if lprint
            if lprint:
                _write_refine(ns, seglist, ali1, s_beg, s_end,
                              segmentrange, idom, est, string,
                              outputunit_lines)

            return flag, est, ali, ali1


def _write_refine(ns, seglist, ali1, s_beg, s_end, segmentrange,
                  idom, est, string, outputunit_lines):
    """Write refine output line."""
    a1 = np.zeros(ns, dtype=np.int32)
    a2 = np.zeros(ns, dtype=np.int32)
    b1 = np.zeros(ns, dtype=np.int32)
    b2 = np.zeros(ns, dtype=np.int32)

    for is_ in range(ns):
        iseg = seglist[is_]
        ires = ali1[iseg]
        if ires != NUL:
            ibeg_key = ires
            if -29 <= ibeg_key < s_beg.shape[0] - 29:
                ibeg = int(s_beg[ibeg_key + 29, iseg])
                iend = int(s_end[ibeg_key + 29, iseg])
            else:
                ibeg = 0
                iend = 0
            a1[is_] = segmentrange[0, iseg] + ibeg
            a2[is_] = segmentrange[1, iseg] - iend
            b1[is_] = ires + ibeg
            b2[is_] = b1[is_] + a2[is_] - a1[is_]
        else:
            a1[is_] = NUL
            a2[is_] = NUL
            b1[is_] = NUL
            b2[is_] = NUL

    # Format: refine<cd1cd2> idom score nseg a1 a2 ... b1 b2 ...
    parts = [f' refine{string}']
    parts.append(f' {idom:>10d}')
    parts.append(f' {est:>10d}')
    parts.append(f' {ns:>10d}')
    for is_ in range(ns):
        parts.append(f' {a1[is_]:>10d}')
        parts.append(f' {a2[is_]:>10d}')
    for is_ in range(ns):
        parts.append(f' {b1[is_]:>10d}')
        parts.append(f' {b2[is_]:>10d}')

    outputunit_lines.append(''.join(parts))


def output1(iseg, nir, ex, cut, ali_nali, ali_start, ali_laststart,
            ali_save, trans, idom, start):
    """Enumerate all singlet candidates above cutoff."""
    seglist_0 = np.array([iseg], dtype=np.int32)
    ali = np.zeros(max(iseg + 1, 1), dtype=np.int32)
    for ir in range(nir):
        iwhere = start[iseg, iseg]
        if iwhere < 0:
            continue
        idx = iwhere + ir
        if 0 <= idx < len(ex) and ex[idx] > cut:
            ali[iseg] = ir
            ali_laststart = saveit(
                idom, ali_nali, ali_start, ali_laststart,
                ali_save, trans, 1, seglist_0, ali)
    return ali_laststart


def output2(aseg, bseg, ex, cut, mi, ci0, idom, ali_nali, ali_start,
            ali_laststart, ali_save, trans, start, lseqtl):
    """Enumerate all doublet combinations above cutoff."""
    iseg = min(aseg, bseg)
    jseg = max(aseg, bseg)
    iistart = start[iseg, iseg]
    jjstart = start[jseg, jseg]
    ijstart = start[iseg, jseg]
    seglist_0 = np.array([aseg, bseg], dtype=np.int32)
    jw = mi[jseg]
    nseg = max(aseg, bseg) + 1
    ali = np.zeros(nseg, dtype=np.int32)

    for ir in range(mi[iseg]):
        ires = ci0[ir, iseg]
        iwhere = ijstart + ires * jw
        for jr in range(mi[jseg]):
            jres = ci0[jr, jseg]
            if lseqtl and trans[jres, jseg] < trans[ires, iseg]:
                continue
            ii_idx = iistart + ires
            jj_idx = jjstart + jres
            ij_idx = iwhere + jres
            if ii_idx >= EXDIM or jj_idx >= EXDIM or ij_idx >= EXDIM:
                continue
            x = ex[ii_idx] + ex[jj_idx]
            if iwhere >= 0:
                x += ex[ij_idx]
            if x > cut:
                ali[iseg] = ires
                ali[jseg] = jres
                ali_laststart = saveit(
                    idom, ali_nali, ali_start, ali_laststart,
                    ali_save, trans, 2, seglist_0, ali)
    return ali_laststart


def saveit(idom, ali_nali, ali_start, ali_laststart, ali_save,
           trans, ns, seglist, ali):
    """Store alignment for parent nodes. Returns updated ali_laststart."""
    # Check degeneracy
    ix = ali_start[idom]
    for i in range(ali_nali[idom]):
        match = True
        for is_ in range(ns):
            iseg = seglist[is_]
            if trans[ali[iseg], iseg] != ali_save[ix + is_ + 1]:
                match = False
                break
        if match:
            return ali_laststart
        ix += ns

    if ali_nali[idom] == 0:
        ali_start[idom] = ali_laststart
    ali_nali[idom] += 1
    for is_ in range(ns):
        iseg = seglist[is_]
        ali_laststart += 1
        if ali_laststart < len(ali_save):
            ali_save[ali_laststart] = trans[ali[iseg], iseg]
    return ali_laststart


def shipup(ni, trans1, idom, ldom, domns, domseglist, mi, ci0, trans,
           idom1, idom2, ali_start, ali_nali, ali_save, bl, nres2, lnul):
    """Build work search space from child node results."""
    nseg = len(mi)
    ni = np.zeros(nseg, dtype=np.int32)
    trans1 = np.zeros_like(trans)

    if idom1 == 0:
        # No children — load all
        _loadall(idom, ni, trans1, mi, ci0, trans, domns, domseglist)
    else:
        if ldom[idom1]:
            _loadali(idom1, domns, domseglist, 0, 0, nres2, ni, trans1,
                     bl, ali_start, ali_nali, ali_save)
        else:
            _loadall(idom1, ni, trans1, mi, ci0, trans, domns, domseglist)

        if ldom[idom2]:
            _loadali(idom2, domns, domseglist, 0, 0, nres2, ni, trans1,
                     bl, ali_start, ali_nali, ali_save)
        else:
            _loadall(idom2, ni, trans1, mi, ci0, trans, domns, domseglist)

    # Add NULs if missing
    if lnul:
        for is_ in range(domns[idom]):
            iseg = domseglist[is_, idom] - 1  # 0-based
            has_nul = False
            for ir in range(ni[iseg]):
                if trans1[ir, iseg] == NUL:
                    has_nul = True
                    break
            if not has_nul:
                trans1[ni[iseg], iseg] = NUL
                ni[iseg] += 1

    return ni, trans1


def _loadall(idom, ni, trans1, mi, ci0, trans, domns, domseglist):
    """Load all candidates from mi/trans."""
    for is_ in range(domns[idom]):
        iseg = domseglist[is_, idom] - 1  # 0-based
        ni[iseg] = mi[iseg]
        for ir in range(ni[iseg]):
            trans1[ir, iseg] = trans[ir, iseg]


def _loadali(idom, domns, domseglist, fuzzN, fuzzC, nres2, ni, trans,
             bl, ali_start, ali_nali, ali_save):
    """Load candidates from saved alignments."""
    stres = -29
    # Build candidate mask
    lcand = {}
    for is_ in range(domns[idom]):
        iseg = domseglist[is_, idom] - 1  # 0-based
        lcand[iseg] = set()

    ix = ali_start[idom]
    for iali in range(ali_nali[idom]):
        for is_ in range(domns[idom]):
            ix += 1
            if ix >= len(ali_save):
                break
            ires = ali_save[ix]
            iseg = domseglist[is_, idom] - 1
            if ires != NUL:
                for i in range(max(stres, ires - fuzzN),
                               min(nres2, ires + fuzzC) + 1):
                    lcand[iseg].add(i)

    for is_ in range(domns[idom]):
        iseg = domseglist[is_, idom] - 1
        ni[iseg] = 0
        for ires in range(stres, nres2 + 1, bl):
            if ires in lcand[iseg]:
                trans[ni[iseg], iseg] = ires
                ni[iseg] += 1


def compress(idom, ni, trans1, mi, ci0, trans, ex, nix, laststart,
             ixstart, domns, domseglist, start):
    """Remap ex[] indices after merging child nodes."""
    nseg = len(mi)
    # Build mapping: old ir -> new ir
    irold = np.zeros((MAXRES0, nseg), dtype=np.int32)
    for is_ in range(domns[idom]):
        iseg = domseglist[is_, idom] - 1  # 0-based
        for ir in range(ni[iseg]):
            ires = 0
            while ires < mi[iseg] and trans[ires, iseg] != trans1[ir, iseg]:
                ires += 1
            if ires >= mi[iseg]:
                return nix, laststart  # error
            irold[ir, iseg] = ires

    # Remap singlets and doublets
    for is_ in range(domns[idom]):
        iseg = domseglist[is_, idom] - 1
        iwhere = start[iseg, iseg]
        if iwhere >= 0:
            for ir in range(min(EXDIM - iwhere, ni[iseg])):
                old_idx = iwhere + irold[ir, iseg]
                new_idx = iwhere + ir
                if 0 <= old_idx < EXDIM and 0 <= new_idx < EXDIM:
                    ex[new_idx] = ex[old_idx]

        for js in range(is_):
            jseg = domseglist[js, idom] - 1
            aseg = min(iseg, jseg)
            bseg = max(iseg, jseg)
            if start[aseg, bseg] >= 0:
                for ir in range(ni[aseg]):
                    iwhere_new = start[aseg, bseg] + ir * ni[bseg]
                    iwhere_old = start[aseg, bseg] + irold[ir, aseg] * mi[bseg]
                    for jr in range(ni[bseg]):
                        old_idx = iwhere_old + irold[jr, bseg]
                        new_idx = iwhere_new + jr
                        if 0 <= old_idx < EXDIM and 0 <= new_idx < EXDIM:
                            ex[new_idx] = ex[old_idx]

    # Update mi, ci0, trans
    for is_ in range(domns[idom]):
        iseg = domseglist[is_, idom] - 1
        for ir in range(ni[iseg]):
            trans[ir, iseg] = trans1[ir, iseg]
            ci0[ir, iseg] = ir
        mi[iseg] = ni[iseg]

    # Left-shift ex chunks (0-based: data at positions start..start+count-1)
    new_laststart = 0
    for ix in range(nix):
        iseg = ixstart[0, ix]
        jseg = ixstart[1, ix]
        old_start = start[iseg, jseg]
        start[iseg, jseg] = new_laststart
        start[jseg, iseg] = new_laststart
        if iseg == jseg:
            count = mi[iseg]
        else:
            count = mi[iseg] * mi[jseg]
        for i in range(count):
            src = old_start + i
            dst = new_laststart + i
            if 0 <= src < EXDIM and 0 <= dst < EXDIM:
                ex[dst] = ex[src]
        new_laststart += count

    return nix, new_laststart


def loadstack(pstack, ali_start, ali_nali, ali_save,
              mi, ci0, trans, idom, domns, domseglist, ex,
              start, lfix, idom1, idom2, ns, seglist, cutoff,
              lseqtl, lnul):
    """Initialize priority queue from child alignments."""
    nseg = len(mi)
    _debug = False  # TEMP

    # Get fixed segments from children
    slist = []
    tlist = []
    if idom1 > 0:
        for is_ in range(domns[idom1]):
            iseg_1 = domseglist[is_, idom1]
            if lfix[iseg_1 - 1, idom1]:
                slist.append(is_)
    if idom2 > 0:
        for is_ in range(domns[idom2]):
            iseg_1 = domseglist[is_, idom2]
            if lfix[iseg_1 - 1, idom2]:
                tlist.append(is_)

    nali1 = len(slist)
    nali2 = len(tlist)
    ia1 = max(1, ali_nali[idom1] if nali1 > 1 else 1)
    ia2 = max(1, ali_nali[idom2] if nali2 > 1 else 1)

    if _debug:
        print(f'    [loadstack] idom={idom}, idom1={idom1}, idom2={idom2}')
        print(f'    [loadstack] slist={slist}, tlist={tlist}, nali1={nali1}, nali2={nali2}')
        print(f'    [loadstack] ia1={ia1}, ia2={ia2}, cutoff={cutoff}')
        print(f'    [loadstack] ali_nali[{idom1}]={ali_nali[idom1]}, ali_nali[{idom2}]={ali_nali[idom2]}')

    # Build initial search space
    ni_work = np.zeros(nseg, dtype=np.int32)
    ci_work = np.zeros((MAXRES0, nseg), dtype=np.int32)
    for is_ in range(ns):
        iseg = seglist[is_]
        ni_work[iseg] = mi[iseg]
        for ir in range(mi[iseg]):
            ci_work[ir, iseg] = ci0[ir, iseg]

    if _debug:
        print(f'    [loadstack] ni_work={[ni_work[seglist[s]] for s in range(ns)]}')
        print(f'    [loadstack] segs={[seglist[s] for s in range(ns)]}')

    for iali in range(ia1):
        # Set fixed segments from child 1
        if nali1 > 1:
            for k in range(nali1):
                iseg_1 = domseglist[slist[k], idom1]
                iseg = iseg_1 - 1
                ni_work[iseg] = 0
                l = ali_save[ali_start[idom1] + iali * domns[idom1] + slist[k] + 1]
                if l != NUL:
                    _selectcand(l, iseg, trans, ni_work, ci_work, mi)
                if lnul:
                    _selectcand(NUL, iseg, trans, ni_work, ci_work, mi)

        for jali in range(ia2):
            # Set fixed segments from child 2
            if nali2 > 1:
                for k in range(nali2):
                    iseg_1 = domseglist[tlist[k], idom2]
                    iseg = iseg_1 - 1
                    ni_work[iseg] = 0
                    l = ali_save[ali_start[idom2] + jali * domns[idom2] + tlist[k] + 1]
                    if l != NUL:
                        _selectcand(l, iseg, trans, ni_work, ci_work, mi)
                    if lnul:
                        _selectcand(NUL, iseg, trans, ni_work, ci_work, mi)

            # Push state
            ess, est = get_ess(ns, seglist, ni_work, ci_work, ex, start,
                               mi, trans, lseqtl)
            if _debug:
                print(f'    [loadstack] iali={iali}, jali={jali}: est={est}, cutoff={cutoff}, push={est > cutoff}')
                print(f'    [loadstack] ni_work after fixed: {[ni_work[seglist[s]] for s in range(ns)]}')
            if est > cutoff:
                candidates = {}
                for is_ in range(ns):
                    seg = seglist[is_]
                    candidates[seg] = [ci_work[k, seg]
                                       for k in range(ni_work[seg])]
                pstack.push(SearchState(est, candidates))
                if _debug:
                    print(f'    [loadstack] pushed state, pstack size={pstack.size}')


def _selectcand(l, iseg, trans, ni, ci, mi):
    """Select candidate l for segment iseg."""
    for ir in range(mi[iseg]):
        if trans[ir, iseg] == l:
            ni[iseg] += 1
            ci[ni[iseg] - 1, iseg] = ir
            return


def findhomolog(idom1, ali1, ali_nali1, ali_start1, ali_save1,
                ali2, domns, domseglist, lfix1):
    """Find matching segments in child alignments."""
    if idom1 <= 0:
        return

    for is_ in range(domns[idom1]):
        iseg = domseglist[is_, idom1] - 1  # 0-based
        ali2[iseg] = NUL

    for iali in range(ali_nali1[idom1]):
        ix = ali_start1[idom1] + iali * domns[idom1]
        found = True
        for is_ in range(domns[idom1]):
            iseg = domseglist[is_, idom1] - 1
            if lfix1[iseg, idom1]:
                x = ali_save1[ix + is_ + 1]
                if x != NUL:
                    if x < ali1[iseg] or x > ali1[iseg] + 9:
                        found = False
                        break
                ali2[iseg] = x
        if found:
            return


def perresiduescore(ns, seglist, cutoff, lnul, ali1, upper, lower,
                    segmentrange, dist, nres1, dist2, nres2, dist2sum,
                    ss, nseg, idom, domns, domseglist, lseqtl, string,
                    ali2, dist1sum, minseglen, segmentrange0,
                    outputunit_lines, scoretable):
    """Refine alignment at single-residue resolution.

    Takes 10-block alignment in ali1 and refines to individual residue positions.
    """
    # Initialize fresh scoring state
    lexdone, start = initlexdonestart(nseg, ss)
    ex = np.zeros(EXDIM, dtype=np.int32)
    laststart = 0
    nix = 0
    ixstart = np.zeros((2, nseg * nseg), dtype=np.int32)
    trans = np.zeros((MAXRES0, nseg), dtype=np.int32)
    ni = np.zeros(nseg, dtype=np.int32)
    ci = np.zeros((MAXRES0, nseg), dtype=np.int32)
    mi_local = np.zeros(nseg, dtype=np.int32)
    s_beg = np.zeros((nres2 + 30, nseg), dtype=np.int16)
    s_end = np.zeros((nres2 + 30, nseg), dtype=np.int16)

    # Build per-residue search space: ±9 residues around 10-block position
    for is_ in range(ns):
        iseg = seglist[is_]
        i = 0
        if ali1[iseg] != NUL:
            for ires in range(ali1[iseg], min(ali1[iseg] + 10, nres2 + 1)):
                ci[i, iseg] = i
                trans[i, iseg] = ires
                i += 1
        if lnul:
            ci[i, iseg] = i
            trans[i, iseg] = NUL
            i += 1
        ni[iseg] = i
        mi_local[iseg] = i

    # Overwrite fixed segments with homolog alignment
    for is_ in range(ns):
        iseg = seglist[is_]
        if ali2[iseg] != NUL:
            i = 0
            ci[i, iseg] = i
            trans[i, iseg] = ali2[iseg]
            i += 1
            if lnul:
                ci[i, iseg] = i
                trans[i, iseg] = NUL
                i += 1
            ni[iseg] = i
            mi_local[iseg] = i

    # Compute scores with bl=1 (single residue resolution)
    nix, laststart = update_ex(
        mi_local, trans, upper, lower, segmentrange, idom, domns,
        domseglist, ex, dist, nres1, dist2, nres2, dist2sum,
        nix, laststart, ixstart, 1, start, lexdone, lseqtl,
        s_beg, s_end, dist1sum, minseglen, scoretable)

    # Build priority stack — separate first segment candidates
    pstack = PriorityStack()
    ess, est_total = get_ess(ns, seglist, ni, ci, ex, start,
                             mi_local, trans, lseqtl)
    firstest = float(max(est_total, 1))

    iseg_first = seglist[0]
    for i in range(mi_local[iseg_first]):
        ni[iseg_first] = 1
        ci[0, iseg_first] = i
        ess_i, est_i = get_ess(ns, seglist, ni, ci, ex, start,
                                mi_local, trans, lseqtl)
        if est_i > cutoff:
            candidates = {}
            for is_ in range(ns):
                seg = seglist[is_]
                candidates[seg] = [ci[k, seg] for k in range(ni[seg])]
            pstack.push(SearchState(est_i, candidates))
    # Restore ni for first seg
    ni[iseg_first] = mi_local[iseg_first]

    # Initialize ali1 to NUL
    for is_ in range(ns):
        ali1[seglist[is_]] = NUL

    # Run dostack with lprint=True
    flag, est, ali_out, ali1_out = dostack(
        pstack, mi_local, ns, seglist, cutoff, ex,
        trans, start, idom, lseqtl, s_beg, s_end,
        segmentrange0, string, outputunit_lines)

    return est, ali1_out
