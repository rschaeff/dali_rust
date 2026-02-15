"""Scoring functions for PARSI — weights, score table, segment scores."""

import numpy as np
from .parsizes import NUL, INFINIT, EXDIM, MAXRES0, MAXSEG


def _nint(x):
    """Fortran-compatible nint: round half away from zero."""
    return int(np.trunc(x + np.copysign(0.5, x)))


def weights():
    """Compute Gaussian envelope weights.

    weight[i] = nint(100*exp(-(i/20)^2)), clipped to 0 if <5.
    Returns array of 1001 ints (index 0 unused, 1..1000).
    Uses single precision to match Fortran.
    """
    w = np.zeros(1001, dtype=np.int32)
    enveloperadius = np.float32(20.0)
    x = np.float32(1.0) / (enveloperadius * enveloperadius)
    for i in range(1, 1001):
        ii = np.float32(i)
        w[i] = _nint(float(np.float32(100) * np.exp(
            np.float32(-x / np.float32(100) * ii * ii))))
        if w[i] < 5:
            w[i] = 0
    return w


def fillscoretable(weight):
    """Compute 160x400 score lookup table.

    scoretable[b, a] = nint(100*weight[a]*(0.20 - abs(x-y)/x))
    where a indexes dist in protein 1 (a/10.0 Å), b indexes binned dist in protein 2.

    Binning: b<100 → y=b*0.1; 100≤b<125 → y=10+(b-100)*0.4; b≥125 → y=b-125+20

    Uses single precision to match Fortran.
    Returns int32 array of shape (161, 401) — 1-based indexing.
    """
    st = np.zeros((161, 401), dtype=np.int32)
    for b in range(1, 161):
        if b < 100:
            y = np.float32(b) * np.float32(0.1)
        elif b < 125:
            y = np.float32(10.0) + np.float32(b - 100) * np.float32(0.4)
        else:
            y = np.float32(b - 125 + 20.0)
        for a in range(1, 401):
            x = np.float32(a) / np.float32(10.0)
            val = np.float32(100.0) * np.float32(weight[a]) * (
                np.float32(0.20) - np.float32(abs(x - y)) / x)
            st[b, a] = _nint(float(val))
    return st


def selfscore(nseg, segmentrange, dist, weight):
    """Compute segment-segment self-scores using Gaussian envelope.

    ss[i,j] = sum over (ires in seg_i, jres in seg_j) of weight[dist[ires,jres]] * 20
    Cutoff: if ss < 42*lali^2 then ss = 0.

    Returns (nseg, nseg) int32 array, 0-based segment indexing.
    """
    ss = np.zeros((nseg, nseg), dtype=np.int32)
    for iseg in range(nseg):
        for jseg in range(iseg + 1):
            s = 0
            for ires in range(segmentrange[0, iseg] - 1, segmentrange[1, iseg]):
                for jres in range(segmentrange[0, jseg] - 1, segmentrange[1, jseg]):
                    p = dist[ires, jres]
                    s += weight[p] * 20
            lali = ((segmentrange[1, iseg] - segmentrange[0, iseg] + 1) +
                    (segmentrange[1, jseg] - segmentrange[0, jseg] + 1))
            if s < 42 * lali * lali:
                s = 0
            ss[iseg, jseg] = s
            ss[jseg, iseg] = s
    return ss


def getupperlower(dist, nseg, segmentrange):
    """Compute 70%/130% distance-sum bounds per segment pair.

    Returns lower, upper arrays of shape (nseg, nseg).
    """
    lower = np.zeros((nseg, nseg), dtype=np.int32)
    upper = np.zeros((nseg, nseg), dtype=np.int32)
    for iseg in range(nseg):
        for jseg in range(iseg + 1):
            s = 0
            for ires in range(segmentrange[0, iseg] - 1, segmentrange[1, iseg]):
                for jres in range(segmentrange[0, jseg] - 1, segmentrange[1, jseg]):
                    s += dist[ires, jres]
            lower[iseg, jseg] = 70 * s
            lower[jseg, iseg] = lower[iseg, jseg]
            upper[iseg, jseg] = 130 * s
            upper[jseg, iseg] = upper[iseg, jseg]
    return lower, upper


def setcut(ndom, domns, domseglist, segmentrange):
    """Compute score cutoffs per domain node.

    cut[idom] = 9000*(polynomial in lali), capped at lali=200.
    """
    cut = np.zeros(ndom + 1, dtype=np.int32)
    for idom in range(1, ndom + 1):
        ns = domns[idom]
        lali = 0
        for i in range(ns):
            seg = domseglist[i, idom]
            lali += segmentrange[1, seg - 1] - segmentrange[0, seg - 1] + 1
        lali = min(lali, 200)
        l = np.float32(lali)
        x = np.float32(9000.0) * (np.float32(0.83259) + np.float32(0.11186) * l
                       + np.float32(3.3537e-5) * l * l * l
                       + np.float32(1.475e-3) * l * l
                       - np.float32(1.579e-7) * l * l * l * l)
        cut[idom] = max(int(x), 0)
    return cut


def flex(ss, nseg, ndom, domns, domseglist):
    """Determine fixed/free segments per domain node (density-based).

    Returns lfix, lfix1 boolean arrays of shape (nseg, ndom+1).
    lfix uses cutoff 0.70, lfix1 uses cutoff 0.90.
    """
    lfix = np.zeros((nseg, ndom + 1), dtype=bool)
    lfix1 = np.zeros((nseg, ndom + 1), dtype=bool)

    # Compute si and sf
    si = np.zeros(nseg, dtype=np.int32)
    sf = np.zeros((nseg, nseg), dtype=np.float32)
    for iseg in range(nseg):
        for jseg in range(nseg):
            if jseg != iseg:
                si[iseg] += ss[iseg, jseg]
        if si[iseg] == 0:
            si[iseg] = 1
        for jseg in range(nseg):
            if jseg != iseg:
                sf[iseg, jseg] = np.float32(ss[iseg, jseg]) / np.float32(si[iseg])

    for idom in range(ndom, 0, -1):
        ns = domns[idom]
        seglist = [domseglist[i, idom] for i in range(ns)]  # 1-based seg indices
        # getlfix with cutoff 0.70
        _getlfix(idom, lfix, 0.70, ns, seglist, sf)
        # getlfix with cutoff 0.90
        _getlfix(idom, lfix1, 0.90, ns, seglist, sf)
        # If domain size > 2 and no segment fixed, fix first
        if ns > 2:
            n = sum(1 for seg in seglist if lfix[seg - 1, idom])
            if n == 0:
                lfix[seglist[0] - 1, idom] = True

    return lfix, lfix1


def _getlfix(idom, lfix, fixcutoff, ns, seglist, sf):
    """Mark segments as fixed if their density score >= fixcutoff."""
    for is_ in range(ns):
        iseg = seglist[is_] - 1  # 0-based
        x = np.float32(0.0)
        for js in range(ns):
            if js != is_:
                jseg_1based = seglist[js]
                x = np.float32(x + sf[iseg, jseg_1based - 1])
        lfix[iseg, idom] = (x >= np.float32(fixcutoff))


def setminseglen(secstr, nseg):
    """Min segment length: 6 for E, 8 for H.

    Returns array of nseg ints (0-based).
    """
    from .parsizes import LENE, LENH
    minseglen = np.zeros(nseg, dtype=np.int32)
    for iseg in range(nseg):
        minseglen[iseg] = LENE
        if secstr[iseg] == 'H':
            minseglen[iseg] = LENH
    return minseglen


def setngap(segmentrange, minseglen, nseg):
    """Compute N-terminal gap positions per segment.

    ngap[iseg] = (lres + bl - 1) / bl where lres = (a2-a1+1-minlen), clamped to >= -29.
    """
    from .parsizes import BL
    ngap = np.zeros(nseg, dtype=np.int32)
    for iseg in range(nseg):
        a2 = segmentrange[1, iseg]
        a1 = segmentrange[0, iseg]
        minlen = minseglen[iseg]
        lres = a2 - a1 + 1 - minlen
        if lres < -29:
            lres = -29
        ngap[iseg] = (lres + BL - 1) // BL
    return ngap


def initlexdonestart(nseg, ss):
    """Initialize lexdone and start arrays.

    lexdone[i,j] = True if ss[i,j] == 0 (nothing to compute).
    start[i,j] = -INFINIT (unassigned).
    """
    lexdone = np.zeros((nseg, nseg), dtype=bool)
    start = np.full((nseg, nseg), -INFINIT, dtype=np.int32)
    for iseg in range(nseg):
        for jseg in range(nseg):
            if ss[iseg, jseg] <= 0:
                lexdone[iseg, jseg] = True
    return lexdone, start


def trimtable(table, ibeg, iend, seglen, minlen):
    """Trim destabilizing end residues from score table.

    Iteratively removes N/C-terminal residues with most negative row sum
    until no negative end row sums remain or min segment length reached.

    Returns (ibeg, iend, totscore).
    """
    while True:
        if seglen - ibeg - iend <= minlen:
            break
        # Compute row sums
        rowsum = np.zeros(seglen + 1, dtype=np.int32)
        for i in range(1 + ibeg, seglen - iend + 1):
            for j in range(1 + ibeg, seglen - iend + 1):
                rowsum[i] += table[i, j]
        # Remove most negative end
        if rowsum[1 + ibeg] < rowsum[seglen - iend]:
            if rowsum[1 + ibeg] < 0:
                ibeg += 1
            else:
                break
        else:
            if rowsum[seglen - iend] < 0:
                iend += 1
            else:
                break

    # Compute total score
    totscore = 0
    for i in range(1 + ibeg, seglen - iend + 1):
        for j in range(1 + ibeg, seglen - iend + 1):
            totscore += table[i, j]

    return ibeg, iend, totscore


def singletex(iseg, upp, low, a1, a2, nir, trans, dist2sum, nres2,
              dist2, nres1, dist, ex, bl, start, s_beg, s_end, dist1sum,
              minlen, scoretable):
    """Compute singlet scores for all candidates of a segment.

    For each candidate ir: try all bl shift positions, compute per-residue-pair
    score using scoretable, trim ends, remember best score.
    """
    iwhere = start[iseg, iseg]
    if iwhere < 0:
        return

    seglen = a2 - a1 + 1
    if seglen > 100:
        return

    for ir in range(nir):
        transires = trans[ir, iseg]

        if transires == NUL:
            x = 0
        else:
            x = -INFINIT
            for l in range(bl):
                ibeg = 0
                iend = 0
                # Beyond C-terminus
                if transires + l > nres2 - minlen + 1:
                    t = -INFINIT
                    # Still record s_beg/s_end
                    if (transires != NUL and transires + l <= nres2
                            and transires + l >= -29):
                        s_beg[transires + l + 29, iseg] = ibeg
                        s_end[transires + l + 29, iseg] = iend
                    continue

                # N- or C-terminal gap
                if transires + l < 1:
                    ibeg = 1 - (transires + l)
                if transires + l + seglen - 1 > nres2:
                    iend = seglen - (nres2 - (transires + l) + 1)

                if ibeg > seglen - minlen:
                    t = -INFINIT
                    ibeg_save = INFINIT
                    iend_save = INFINIT
                else:
                    # Calculate table of per-residue scores
                    table = np.zeros((101, 101), dtype=np.int32)
                    for i in range(ibeg, seglen - iend):
                        for j in range(ibeg, seglen - iend):
                            # Fortran: p=(a1+i-1)*nres1+a1+j, using 1-based flat index
                            # In our 0-based arrays:
                            p_row = (a1 - 1) + i  # 0-based row in dist
                            p_col = (a1 - 1) + j  # 0-based col in dist
                            q_row = (transires + l - 1) + i  # 0-based row in dist2
                            q_col = (transires + l - 1) + j  # 0-based col in dist2
                            if (0 <= q_row < nres2 and 0 <= q_col < nres2
                                    and 0 <= p_row < nres1 and 0 <= p_col < nres1):
                                table[i + 1, j + 1] = scoretable[
                                    dist2[q_row, q_col], dist[p_row, p_col]]

                    # Trim segment ends
                    ibeg, iend, t = trimtable(table, ibeg, iend, seglen, minlen)
                    ibeg_save = ibeg
                    iend_save = iend

                # Remember individual end-gaps
                if (transires != NUL and transires + l <= nres2
                        and transires + l >= -29):
                    s_beg[transires + l + 29, iseg] = ibeg_save if ibeg <= seglen else INFINIT
                    s_end[transires + l + 29, iseg] = iend_save if iend <= seglen else INFINIT

                x = max(t, x)

        # Assign candidate-segment score
        idx = iwhere + ir
        if 0 <= idx < len(ex):
            ex[idx] = x


def segsegscore(iseg, jseg, transires, transjres, a1, a2, b1, b2,
                dist, nres1, dist2, nres2, upp1, low1, dist2sum, bl,
                lseqtl, s_beg, s_end, dist1sum, scoretable):
    """Compute segment-segment doublet score.

    Tries all ishift×jshift combinations (0..bl-1 each).
    Uses distance sum filter and sequential topology enforcement.
    """
    if transires == NUL or transjres == NUL:
        return 0

    s = -INFINIT
    if lseqtl and iseg > jseg and transires < transjres:
        return s
    if lseqtl and iseg < jseg and transires > transjres:
        return s

    lself = (iseg == jseg)
    j1 = transires - a1
    j2 = transjres - b1

    for jshift in range(bl):
        if transjres + jshift > nres2:
            break
        # s_beg/s_end use transjres+jshift as key (can be negative for gaps)
        jbeg_key = transjres + jshift
        if jbeg_key < -29 or jbeg_key > nres2:
            continue
        jbeg = int(s_beg[jbeg_key + 29, jseg])
        jend = int(s_end[jbeg_key + 29, jseg])
        l0 = (jbeg == 0 and jend == 0)
        k2 = j2 + jshift
        d1 = b1 + k2 - 1 + jbeg  # 1-based
        d2 = b2 + k2 - jend      # 1-based
        e1 = b1 + jbeg - 1       # 1-based
        e2 = b2 - jend            # 1-based

        # Check chain ends
        jfirst = transjres + jshift + jbeg
        jlast = transjres + jshift - jend + b2 - b1
        if jlast > nres2:
            break

        for ishift in range(bl):
            if lself and jshift != ishift:
                continue
            if transires + ishift > nres2:
                continue

            ibeg_key = transires + ishift
            if ibeg_key < -29 or ibeg_key > nres2:
                continue
            ibeg = int(s_beg[ibeg_key + 29, iseg])
            iend = int(s_end[ibeg_key + 29, iseg])
            l1 = (l0 and ibeg == 0 and jbeg == 0)
            k1 = j1 + ishift

            # Check chain ends
            ifirst = transires + ishift + ibeg
            ilast = transires + ishift - iend + a2 - a1
            if ilast > nres2:
                continue

            # Disallow segment overlaps
            if iseg < jseg and ilast >= jfirst:
                continue
            if jseg < iseg and jlast >= ifirst:
                continue

            # Check distance sum
            # Inline computation of ds2
            ds2 = 0
            for c in range(a1 + k1 + ibeg, a2 + k1 - iend + 1):  # 1-based
                if 1 <= c <= nres2:
                    ds2 += (dist2sum[c - 1, d2] - dist2sum[c - 1, d1])

            if not l1:
                ds1 = 0
                for c in range(a1 + ibeg, a2 - iend + 1):  # 1-based
                    if 1 <= c <= nres1:
                        ds1 += (dist1sum[c - 1, e2] - dist1sum[c - 1, e1])
                low = int(np.float32(0.7) * np.float32(ds1))
                upp = int(np.float32(1.3) * np.float32(ds1))
                if ds2 <= low or ds2 >= upp:
                    continue
            else:
                if ds2 <= low1 or ds2 >= upp1:
                    continue

            # Precise calculation
            x = 0
            for a in range(a1 + ibeg, a2 - iend + 1):  # 1-based
                p_row = a - 1  # 0-based
                q_row = a + k1 - 1  # 0-based
                for b in range(b1 + jbeg, b2 - jend + 1):  # 1-based
                    b_col = b - 1  # 0-based
                    q_col = b + k2 - 1 + 1 - 1  # = b + k2 - 1, 0-based
                    # Actually: Fortran q+b where q = (a1+ibeg+k1-1)*nres2+k2
                    # and b ranges from b1+jbeg to b2-jend
                    # so dist2 index = (q_row, transjres+jshift + (b-b1) + jbeg - 1 + (b-b1-jbeg)...)
                    # Let me recalculate from Fortran:
                    # p = (a1+ibeg-1)*nres1, then p+b = that + b (1-based flat)
                    # So dist[p_row, b_col] where p_row = a-1, b_col = b-1
                    # q = (a1+ibeg+k1-1)*nres2+k2, then q+b
                    # q_row_fortran = a1+ibeg+k1-1 = a+k1-1 (after loop increment)
                    # Wait, the Fortran loop: p starts at (a1+ibeg-1)*nres1, increments by nres1
                    # q starts at (a1+ibeg+k1-1)*nres2+k2, increments by nres2
                    # For each a (1-based from a1+ibeg to a2-iend):
                    #   dist row in protein1: a-1 (0-based)
                    #   dist2 row: a + k1 - 1 (0-based, since k1 = j1 + ishift = transires-a1+ishift)
                    #     = a + transires - a1 + ishift - 1 = transires + ishift + (a - a1) - 1
                    #     0-based: transires + ishift - 1 + (a - a1)
                    #   dist2 col for b: k2 + b = j2 + jshift + b = transjres - b1 + jshift + b
                    #     0-based: transjres + jshift - 1 + (b - b1)
                    d2_row = transires + ishift - 1 + (a - a1)  # 0-based
                    d2_col = transjres + jshift - 1 + (b - b1)  # 0-based
                    if (0 <= d2_row < nres2 and 0 <= d2_col < nres2
                            and dist2[d2_row, d2_col] >= 1
                            and dist[p_row, b_col] >= 1):
                        x += scoretable[dist2[d2_row, d2_col], dist[p_row, b_col]]
            if x > s:
                s = x

    return s


def doubletex(aseg, bseg, nir, njr, upp, low, a1, a2, b1, b2,
              trans, dist2sum, nres2, dist2, nres1, dist, ex, bl,
              start, lseqtl, s_beg, s_end, dist1sum, scoretable):
    """Compute doublet scores for all candidate pairs of two segments."""
    ijstart = start[aseg, bseg]
    if ijstart < 0:
        return

    for ir in range(nir):
        iwhere = ijstart + ir * njr
        transires = trans[ir, aseg]
        if transires == NUL:
            for jr in range(njr):
                idx = iwhere + jr
                if 0 <= idx < len(ex):
                    ex[idx] = 0
        else:
            for jr in range(njr):
                transjres = trans[jr, bseg]
                if transjres == NUL:
                    x = 0
                else:
                    x = segsegscore(
                        aseg, bseg, transires, transjres,
                        a1, a2, b1, b2, dist, nres1, dist2, nres2,
                        upp, low, dist2sum, bl, lseqtl, s_beg, s_end,
                        dist1sum, scoretable)
                idx = iwhere + jr
                if 0 <= idx < len(ex):
                    ex[idx] = x


def update_ex(mi, trans, upper, lower, segmentrange, idom, domns,
              domseglist, ex, dist, nres1, dist2, nres2, dist2sum,
              nix, laststart, ixstart, bl, start, lexdone, lseqtl,
              s_beg, s_end, dist1sum, minseglen, scoretable):
    """Compute scores for new segment pairs in domain idom.

    Updates ex[], start[], lexdone[] in-place. Returns updated (nix, laststart).
    """
    ns = domns[idom]
    for is_ in range(ns):
        iseg = domseglist[is_, idom] - 1  # 0-based
        if not lexdone[iseg, iseg]:
            nix, laststart = _setstart1(iseg, iseg, nix, laststart, ixstart,
                                        mi, start, lexdone)
            singletex(iseg, upper[iseg, iseg], lower[iseg, iseg],
                      segmentrange[0, iseg], segmentrange[1, iseg],
                      mi[iseg], trans, dist2sum, nres2,
                      dist2, nres1, dist, ex, bl, start, s_beg, s_end,
                      dist1sum, minseglen[iseg], scoretable)

        for js in range(is_):
            jseg = domseglist[js, idom] - 1  # 0-based
            if not lexdone[iseg, jseg]:
                aseg = min(iseg, jseg)
                bseg = max(iseg, jseg)
                nix, laststart = _setstart1(aseg, bseg, nix, laststart,
                                            ixstart, mi, start, lexdone)
                doubletex(aseg, bseg, mi[aseg], mi[bseg],
                          upper[aseg, bseg], lower[aseg, bseg],
                          segmentrange[0, aseg], segmentrange[1, aseg],
                          segmentrange[0, bseg], segmentrange[1, bseg],
                          trans, dist2sum, nres2, dist2, nres1, dist,
                          ex, bl, start, lseqtl, s_beg, s_end, dist1sum,
                          scoretable)

    return nix, laststart


def _setstart1(iseg, jseg, nix, laststart, ixstart, mi, start, lexdone):
    """Set start position for segment pair in ex[] array."""
    ixstart[0, nix] = iseg
    ixstart[1, nix] = jseg
    start[iseg, jseg] = laststart
    start[jseg, iseg] = laststart
    if iseg == jseg:
        laststart = laststart + mi[iseg]
    else:
        laststart = laststart + mi[iseg] * mi[jseg]
    if laststart > EXDIM:
        laststart = 0
    lexdone[iseg, jseg] = True
    lexdone[jseg, iseg] = True
    nix += 1
    return nix, laststart


def get_ess(ns, seglist, ni, ci, ex, start, mi, trans, lseqtl):
    """Compute segment score contributions (ess) and total estimate.

    Returns (ess, est) where ess is (nseg, nseg) array and est is total estimate.
    """
    nseg = len(mi)
    ess = np.zeros((nseg, nseg), dtype=np.int32)
    est = 0
    for is_ in range(ns):
        iseg = seglist[is_]
        e = get_estimate(iseg, iseg, ni, ci, ex, start, mi, trans, lseqtl)
        ess[iseg, iseg] = e
        est += e
        for js in range(is_):
            jseg = seglist[js]
            e = get_estimate(iseg, jseg, ni, ci, ex, start, mi, trans, lseqtl)
            e = e + e  # cross-terms added twice
            ess[iseg, jseg] = e
            ess[jseg, iseg] = e
            est += e
    return ess, est


def get_estimate(iseg, jseg, ni, ci, ex, start, mi, trans, lseqtl):
    """Get maximum score estimate for a segment pair."""
    if iseg < jseg:
        aseg, bseg = iseg, jseg
    else:
        aseg, bseg = jseg, iseg

    if start[iseg, jseg] < 0:
        est = 0
        if lseqtl and aseg != bseg and ni[aseg] >= 1 and ni[bseg] >= 1:
            # Check if sequential topology violated
            ares = trans[ci[ni[aseg] - 1, aseg], aseg]
            if ares == NUL:
                return 0
            bres = trans[ci[ni[bseg] - 1, bseg], bseg]
            if bres == NUL:
                return 0
            ares = trans[ci[0, aseg], aseg]
            if ares > bres:
                est = -INFINIT
        return est
    else:
        est = -INFINIT
        iwhere0 = start[aseg, bseg]
        if aseg == bseg:
            for i in range(ni[aseg]):
                idx = iwhere0 + ci[i, aseg]
                if idx >= len(ex):
                    break
                x = ex[idx]
                if x > est:
                    est = x
        else:
            jw = mi[bseg]
            for i in range(ni[aseg]):
                iwhere = iwhere0 + ci[i, aseg] * jw
                for j in range(ni[bseg]):
                    idx = iwhere + ci[j, bseg]
                    if idx >= len(ex):
                        break
                    x = ex[idx]
                    if x > est:
                        est = x
        return est
