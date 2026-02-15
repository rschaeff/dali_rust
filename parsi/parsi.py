"""PARSI orchestration — dowork_parsi and run_parsi entry points."""

import numpy as np
from pathlib import Path

from .parsizes import NUL, BL, STARTSIZE, LENE, LENH
from .protein_io import (parsireadproteindata, treehack, compressca,
                          getdist, hackdist, getdist2, getdist2sum)
from .scoring import (weights, fillscoretable, selfscore, getupperlower,
                       setcut, flex, setminseglen, setngap,
                       initlexdonestart)
from .domain_tree import setldom, init_searchspace, init_ci_ni
from .align import align


# Module-level cached scoring tables
_weight = None
_scoretable = None


def init_parsi():
    """Initialize PARSI scoring tables (called once)."""
    global _weight, _scoretable
    _weight = weights()
    _scoretable = fillscoretable(_weight)


def dowork_parsi(cd1, cd2, dat_dir1, dat_dir2=None, lfirstonly=True,
                 cd1_cache=None):
    """Process a single protein pair through PARSI.

    Args:
        cd1: Query protein code (e.g. '101mA')
        cd2: Target protein code (e.g. '1a00A')
        dat_dir1: Path to .dat files for protein 1
        dat_dir2: Path to .dat files for protein 2 (defaults to dat_dir1)
        lfirstonly: If True, limit to 5 alignments per domain
        cd1_cache: Dict to cache protein 1 data between calls

    Returns:
        List of refine output lines.
    """
    global _weight, _scoretable
    if _weight is None:
        init_parsi()

    if dat_dir2 is None:
        dat_dir2 = dat_dir1

    output_lines = []

    # === Phase 1: Load protein 1 (cached if same cd1) ===
    if cd1_cache is not None and cd1_cache.get('code') == cd1:
        p1 = cd1_cache['prot']
        dist1 = cd1_cache['dist']
        ss = cd1_cache['ss']
        upper = cd1_cache['upper']
        lower = cd1_cache['lower']
        lfix = cd1_cache['lfix']
        lfix1 = cd1_cache['lfix1']
        cut = cd1_cache['cut']
        minseglen = cd1_cache['minseglen']
        ngap = cd1_cache['ngap']
        dist1sum = cd1_cache['dist1sum']
        ldom = cd1_cache['ldom']
    else:
        dat_path = Path(dat_dir1) / f'{cd1}.dat'
        p1 = parsireadproteindata(dat_path)

        # Skip proteins with too few SSEs
        if p1.nseg <= 2 or p1.ndom < p1.nseg:
            return output_lines
        if p1.na * p1.na + p1.nb * p1.nb == 0:
            return output_lines

        treehack(p1)

        # Save original ranges, then compress
        compressca(p1)
        if p1.nres == 0:
            return output_lines

        # Distance matrix
        dist1 = getdist(p1)

        # Self-scores
        ss = selfscore(p1.nseg, p1.segmentrange, dist1, _weight)

        # Upper/lower bounds
        lower, upper = getupperlower(dist1, p1.nseg, p1.segmentrange)

        # Hack distances (cap at 400)
        hackdist(dist1)

        # Flex: fixed/free segments
        lfix, lfix1 = flex(ss, p1.nseg, p1.ndom, p1.domns, p1.domseglist)

        # Score cutoffs
        cut = setcut(p1.ndom, p1.domns, p1.domseglist, p1.segmentrange)

        # Cumulative distance sums for protein 1
        dist1sum = getdist2sum(p1.ca, p1.nres)

        # Domain alignment flags
        ldom = setldom(p1.ndom, p1.domns)

        # Min segment lengths and gap positions
        minseglen = setminseglen(p1.secstr, p1.nseg)
        ngap = setngap(p1.segmentrange, minseglen, p1.nseg)

        # Cache
        if cd1_cache is not None:
            cd1_cache.update({
                'code': cd1,
                'prot': p1,
                'dist': dist1,
                'ss': ss,
                'upper': upper,
                'lower': lower,
                'lfix': lfix,
                'lfix1': lfix1,
                'cut': cut,
                'minseglen': minseglen,
                'ngap': ngap,
                'dist1sum': dist1sum,
                'ldom': ldom,
            })

    # === Phase 2: Load protein 2 ===
    dat_path2 = Path(dat_dir2) / f'{cd2}.dat'
    p2 = parsireadproteindata(dat_path2)

    if p2.nres == 0:
        return output_lines
    if p1.na * p2.na + p1.nb * p2.nb == 0:
        return output_lines

    # Build segment mapping for protein 2
    segment2 = np.zeros(p2.nres, dtype=np.int32)
    for iseg in range(p2.nseg):
        for ires in range(p2.segmentrange[0, iseg] - 1, p2.segmentrange[1, iseg]):
            segment2[ires] = iseg + 1

    # Distance sums and binned distance matrix for protein 2
    dist2sum = getdist2sum(p2.ca, p2.nres)
    dist2 = getdist2(p2.ca, p2.nres)

    # === Phase 3: Initialize search space and run alignment ===
    nseg = p1.nseg

    # Initialize search space
    trans, mi = init_searchspace(nseg, p2.nres, p1.segmentrange, minseglen,
                                 ngap, BL)
    ci0, _ = init_ci_ni(mi, nseg)

    # Build string identifier
    string = f'{cd1:5s}{cd2:5s}'

    # Run alignment
    align(p1, p2.nres, p2.secstr, p2.nseg, segment2,
          dist1, p1.nres, dist2, p2.nres, dist2sum, dist1sum,
          ss, upper, lower, p1.segmentrange, p1.segmentrange0,
          p1.ndom, p1.node_child, p1.domns, p1.domseglist, ldom,
          lfix, lfix1, cut, minseglen, ngap,
          p1.checkrange, p1.checkx,
          mi, ci0, trans, _scoretable, _weight,
          True, True, lfirstonly, string, output_lines)

    return output_lines


def run_parsi(structures, dat_dir, dat_dir2=None):
    """Run PARSI on all pairs of structures.

    Args:
        structures: List of structure codes
        dat_dir: Path to .dat files
        dat_dir2: Path to .dat files for protein 2 (defaults to dat_dir)

    Returns:
        List of all refine output lines.
    """
    init_parsi()
    all_lines = []
    cd1_cache = {}

    for cd1 in structures:
        for cd2 in structures:
            lines = dowork_parsi(cd1, cd2, dat_dir, dat_dir2,
                                 lfirstonly=True, cd1_cache=cd1_cache)
            all_lines.extend(lines)

    return all_lines
