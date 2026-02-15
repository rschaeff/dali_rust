"""Domain tree handling for PARSI."""

import numpy as np
from .parsizes import NUL, BL, MAXRES0, STARTSIZE


def setldom(ndom, domns):
    """Mark domains with >= STARTSIZE segments as 'align'.

    Returns boolean array ldom[0..ndom], with ldom[0]=False.
    """
    # Compute node_size: number of segments (domns already has this)
    ldom = np.zeros(ndom + 1, dtype=bool)
    for idom in range(1, ndom + 1):
        ldom[idom] = (domns[idom] >= STARTSIZE)
    return ldom


def init_searchspace(nseg, nres2, segmentrange, minseglen, ngap, bl=BL):
    """Initialize search space: enumerate candidate residues for each segment.

    For each segment, candidates are:
    1. Residues in protein 2 at 10-residue steps
    2. N-terminal gap residues (negative indices down to -29)
    3. NUL candidate (unaligned)

    Returns:
        trans: (maxcand, nseg) — maps candidate index to residue number (or NUL)
        mi: (nseg,) — number of candidates per segment
    """
    max_cand = MAXRES0
    trans = np.full((max_cand, nseg), NUL, dtype=np.int32)
    mi = np.zeros(nseg, dtype=np.int32)

    for iseg in range(nseg):
        ir = 0
        # Regular candidates: step through protein 2 at bl intervals
        for ires in range(1, nres2 + 1, bl):
            if ir >= max_cand:
                break
            trans[ir, iseg] = ires
            ir += 1

        # N-terminal gap candidates
        for igap in range(1, ngap[iseg] + 1):
            if ir >= max_cand:
                break
            gapres = 1 - igap * bl
            if gapres < -29:
                gapres = -29
            trans[ir, iseg] = gapres
            ir += 1

        # NUL candidate
        if ir < max_cand:
            trans[ir, iseg] = NUL
            ir += 1

        mi[iseg] = ir

    return trans, mi


def init_ci_ni(mi, nseg):
    """Initialize candidate index arrays — all candidates active.

    ci[ir, iseg] = ir (identity mapping)
    ni[iseg] = mi[iseg]
    """
    max_cand = MAXRES0
    ci = np.zeros((max_cand, nseg), dtype=np.int32)
    ni = np.zeros(nseg, dtype=np.int32)
    for iseg in range(nseg):
        for ir in range(mi[iseg]):
            ci[ir, iseg] = ir
        ni[iseg] = mi[iseg]
    return ci, ni
