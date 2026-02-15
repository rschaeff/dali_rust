"""
FILTER95 module — redundancy filtering with FITZ refinement.

Translated from comparemodules.f dowork_filter95() and filter95fitz.f.

Takes PARSI refine output lines and:
1. Reads domain tree info from .dat files
2. Filters by domain type (keep only '*' and '+')
3. Applies FITZ superposition refinement
4. Scores with distance-based scorefun95
5. Filters by Z-score threshold

Parameters (from serialcompare.f):
    zcut1 = 1.0
    fitzrcut = 4.0
    fitzmaxiter = 3
"""

import numpy as np
from wolf.alignment import getut, filltable_maxsim, nw_maxsim
from wolf.superposition import transrotate


# ---------------------------------------------------------------------------
# FITZ for FILTER95 (takes initial alignment, matching Fortran behavior)
# ---------------------------------------------------------------------------

def fitz95(x, y, ali, rcut, maxiter):
    """
    Iterative superposition refinement — Fortran-compatible version.

    Unlike wolf.alignment.fitz, this takes an initial alignment array and
    checks convergence against it on the first iteration. This matches the
    Fortran fitz() which takes ali as in-out.

    Modifies x in-place.

    Args:
        x: (3, nx) first coordinate set (modified in-place)
        y: (3, ny) second coordinate set
        ali: (nx,) initial alignment, ali[i] = j (0-based) or -1
        rcut: distance cutoff
        maxiter: max iterations

    Returns:
        ali: (nx,) final alignment
        rms: final RMSD
        lali: alignment length
        niter: number of iterations
    """
    nx = x.shape[1]
    ny = y.shape[1]

    if nx < 1 or ny < 1:
        return np.full(max(nx, 0), -1, dtype=np.int32), 0.0, 0, 0

    rms = 0.0
    lali = 0
    niter = 0
    frozen = False

    while niter < maxiter and not frozen:
        niter += 1

        table0 = filltable_maxsim(x, y, rcut)
        ali_old = ali.copy()
        score, ali = nw_maxsim(nx, ny, table0)

        frozen = np.array_equal(ali, ali_old)

        u, t, lali, rms = getut(nx, ali, x, y)
        transrotate(x, u, t)

    return ali, rms, lali, niter


# ---------------------------------------------------------------------------
# .dat reader with domain tree support
# ---------------------------------------------------------------------------

def readproteindata95(dat_path):
    """
    Read .dat file including domain tree info.

    Returns:
        nres: number of residues
        ca: (3, nres) CA coordinates
        ndom: number of domains
        node_type: list[str] of length ndom+1 (1-indexed), node_type[i] for domain i
        node_size: list[int] of length ndom+1 (1-indexed)
        node_nseg: list[int] of length ndom+1 (1-indexed)
        segmentrange: dict {idom: list of (start, end) tuples} (1-based residue indices)
    """
    with open(dat_path) as f:
        lines = f.readlines()

    idx = 0

    # Header: >>>> code nres nseg na nb secstr
    header = lines[idx]
    # Format: (10x,4i5,2x,200a1) — skip 10 chars, read 4 ints (5-wide)
    nres = int(header[10:15])
    nseg = int(header[15:20])
    idx += 1

    # Segment ranges: nseg lines, format (6i10)
    for _ in range(nseg):
        idx += 1

    # CA coordinates: format (10f8.1), 3*nres values
    ca_vals = []
    needed = 3 * nres
    while len(ca_vals) < needed:
        line = lines[idx]
        # Parse 10 floats per line, 8 chars wide
        for k in range(0, len(line.rstrip()), 8):
            chunk = line[k:k+8].strip()
            if chunk:
                ca_vals.append(float(chunk))
                if len(ca_vals) >= needed:
                    break
        idx += 1

    ca = np.zeros((3, nres), dtype=np.float32)
    for i in range(nres):
        ca[0, i] = ca_vals[3*i]
        ca[1, i] = ca_vals[3*i+1]
        ca[2, i] = ca_vals[3*i+2]

    # First domain tree block: >>>> code ndom
    header2 = lines[idx]
    ndom_1 = int(header2[10:15])
    idx += 1

    # Skip first block (segment-based domain tree)
    for _ in range(ndom_1):
        idx += 1

    # Second domain tree block: >>>> code ndom
    header3 = lines[idx]
    ndom = int(header3[10:15])
    idx += 1

    # Read residue-based domain tree
    # Format 530: (i4,1x,a1,1x,4i4,64000i4)
    node_type = [''] * (ndom + 1)   # 1-indexed
    node_size = [0] * (ndom + 1)
    node_nseg = [0] * (ndom + 1)
    segmentrange = {}

    for _ in range(ndom):
        line = lines[idx]
        idx += 1

        # Parse: idom type child1 child2 size nseg ranges...
        idom = int(line[0:4])
        ntype = line[5]
        # Skip child1 (7:11) and child2 (11:15)
        size = int(line[15:19])
        dom_nseg = int(line[19:23])

        node_type[idom] = ntype
        node_size[idom] = size
        node_nseg[idom] = dom_nseg

        # Parse segment ranges: pairs of (start, end)
        ranges = []
        rest = line[23:]
        vals = rest.split()
        for k in range(0, dom_nseg * 2, 2):
            rstart = int(vals[k])
            rend = int(vals[k+1])
            ranges.append((rstart, rend))
        segmentrange[idom] = ranges

    return nres, ca, ndom, node_type, node_size, node_nseg, segmentrange


# ---------------------------------------------------------------------------
# Distance matrix
# ---------------------------------------------------------------------------

def getdist95(ca, nres):
    """
    Compute integer*2 distance matrix (distances * 100).

    Args:
        ca: (3, nres) coordinates
        nres: number of residues

    Returns:
        d: (nres, nres) int16 distance matrix
    """
    d = np.zeros((nres, nres), dtype=np.int16)
    for i in range(nres):
        for j in range(i + 1, nres):
            dx = ca[0, i] - ca[0, j]
            dy = ca[1, i] - ca[1, j]
            dz = ca[2, i] - ca[2, j]
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            val = int(np.rint(100.0 * dist))
            # int16 range: -32768 to 32767
            if val > 32767:
                val = 32767
            d[i, j] = val
            d[j, i] = val
    return d


# ---------------------------------------------------------------------------
# Scoring functions
# ---------------------------------------------------------------------------

def scorefun95(r1, r2):
    """
    Distance-based scoring function.

    Args:
        r1, r2: int16 distance values (distance * 100)

    Returns:
        score: float
    """
    r = (int(r1) + int(r2)) / 200.0
    s = (int(r1) - int(r2)) / 100.0
    if r > 0.01:
        return (0.20 - abs(s) / r) * np.exp(-r * r / 400.0)
    else:
        return 0.20


def gettotscore95(ali, d1, nres1, d2, nres2):
    """
    Compute total score from alignment and distance matrices.

    Args:
        ali: (nres1,) alignment array, ali[i] = j (1-based) or 0
        d1: (nres1, nres1) int16 distance matrix
        d2: (nres2, nres2) int16 distance matrix
        nres1, nres2: sizes

    Returns:
        totscore: float
    """
    # Collect aligned positions
    aligned = []
    for i in range(min(nres1, len(ali))):
        if ali[i] > 0:
            aligned.append(i)

    totscore = 0.0
    for i_idx in range(len(aligned)):
        k = aligned[i_idx]
        for j_idx in range(len(aligned)):
            l = aligned[j_idx]
            q = abs(ali[k]) - 1  # Convert to 0-based for d2
            r = abs(ali[l]) - 1
            if q < 0 or r < 0:
                continue
            if q >= nres2 or r >= nres2:
                continue
            totscore += scorefun95(d1[k, l], d2[q, r])

    return totscore


# ---------------------------------------------------------------------------
# Z-score baseline
# ---------------------------------------------------------------------------

def compute_minscore(node_size_val):
    """
    Compute Z-score baseline mean and sigma for a domain size.

    Returns:
        mean_int: nint(10000 * mean)
        sigma_int: nint(10000 * sigma)
    """
    x = min(float(node_size_val), 400.0)
    mean = 7.9494 + 0.70852 * x + 2.5895e-4 * x * x - 1.9156e-6 * x * x * x
    sigma = max(1.0, 0.50 * mean)
    return int(np.rint(10000 * mean)), int(np.rint(10000 * sigma))


# ---------------------------------------------------------------------------
# Parse PARSI refine lines
# ---------------------------------------------------------------------------

def parse_refine_line(line):
    """
    Parse a PARSI refine output line.

    Format: ' refine<cd1><cd2> idom score nseg a1 a2 ... b1 b2 ...'
    cd1 at positions 7:12, cd2 at 12:17 (0-based Python)

    Returns:
        cd1, cd2, idom, score, nseg, ali (list of ints, nseg*4 values)
    """
    line = line.rstrip()
    if 'refine' not in line:
        return None

    # Find "refine" and extract cd1/cd2
    rpos = line.index('refine')
    cd1 = line[rpos+6:rpos+11]
    cd2 = line[rpos+11:rpos+16]
    rest = line[rpos+16:]

    parts = rest.split()
    idom = int(parts[0])
    score = int(parts[1])
    nseg = int(parts[2])
    ali = [int(x) for x in parts[3:3+nseg*4]]

    return cd1, cd2, idom, score, nseg, ali


# ---------------------------------------------------------------------------
# Main FILTER95
# ---------------------------------------------------------------------------

def dowork_filter95(refine_lines, dat_dir, zcut1=1.0, fitzrcut=4.0,
                    fitzmaxiter=3):
    """
    Process PARSI refine lines through FILTER95.

    Mirrors comparemodules.f dowork_filter95() — the fast-path Z-score
    bypass is COMMENTED OUT in that version, so ALL alignments go through
    FITZ refinement.

    Args:
        refine_lines: list of PARSI refine output lines
        dat_dir: path to .dat files directory (with trailing /)
        zcut1: Z-score cutoff (default 1.0)
        fitzrcut: FITZ distance cutoff (default 4.0)
        fitzmaxiter: FITZ max iterations (default 3)

    Returns:
        output_lines: list of output strings
    """
    output_lines = []

    # Cached state — only cd1/cd2 data persists across calls (matching Fortran)
    # In the module version (comparemodules.f), oldidom is a LOCAL variable
    # initialized to 0 on every call, meaning:
    # 1. The dedup check never triggers (idom is never 0)
    # 2. Domain data (resix/xiser/xca) is always rebuilt
    # 3. xca is always fresh from ca
    old_cd1 = '?????'
    old_cd2 = '?????'

    # cd1 data (cached across lines with same cd1)
    nres1 = 0
    ca1 = None
    d1 = None
    ndom1 = 0
    node_type1 = []
    node_size1 = []
    node_nseg1 = []
    segmentrange1 = {}
    lkeep = {}
    minscore = {}  # {idom: (mean_int, sigma_int)}

    # cd2 data (cached across lines with same cd2)
    nres2 = 0
    ca2 = None
    d2 = None

    for line in refine_lines:
        parsed = parse_refine_line(line)
        if parsed is None:
            continue

        cd1, cd2, idom, score, nseg, ali = parsed

        # Load cd1 data if changed
        if cd1 != old_cd1:
            dat_path = dat_dir + cd1 + '.dat'
            nres1, ca1, ndom1, node_type1, node_size1, node_nseg1, \
                segmentrange1 = readproteindata95(dat_path)
            d1 = getdist95(ca1, nres1)

            lkeep = {}
            for i in range(1, ndom1 + 1):
                lkeep[i] = (node_type1[i] == '*' or node_type1[i] == '+')

            minscore = {}
            for i in range(1, ndom1 + 1):
                minscore[i] = compute_minscore(node_size1[i])

            old_cd1 = cd1

        # Skip non-kept domains
        if not lkeep.get(idom, False):
            continue

        # ALWAYS rebuild domain-local data (oldidom=0 in Fortran module version)
        resix_list = []
        xiser_dict = {}

        for seg_start, seg_end in segmentrange1.get(idom, []):
            for j in range(seg_start, seg_end + 1):
                if j > 0 and j <= nres1:
                    resix_list.append(j)
                    xiser_dict[j] = len(resix_list)  # 1-based

        nx = len(resix_list)
        # Build fresh xca from ca (always rebuilt)
        xca = np.zeros((3, nx), dtype=np.float32)
        for k in range(nx):
            res_idx = resix_list[k] - 1  # 0-based for ca array
            xca[:, k] = ca1[:, res_idx]

        resix = resix_list
        xiser = xiser_dict

        # Load cd2 data if changed
        if cd2 != old_cd2:
            dat_path = dat_dir + cd2 + '.dat'
            nres2, ca2_full, _, _, _, _, _ = readproteindata95(dat_path)
            ca2 = ca2_full
            d2 = getdist95(ca2, nres2)
            old_cd2 = cd2

        # Map PARSI alignment to fitzali in domain-local coords
        # fitzali: 0-based domain-local index -> 1-based cd2 residue (or -1)
        fitzali = np.full(nx, -1, dtype=np.int32)

        for i in range(nseg):
            k = i * 2
            cd1_start = ali[k]      # 1-based, or -99 (NUL)
            cd1_end = ali[k + 1]    # 1-based
            cd2_start = ali[k + nseg * 2]  # 1-based
            # cd2_end = ali[k + nseg * 2 + 1]  # not used directly

            if cd1_start > 0:
                for j in range(cd1_start, cd1_end + 1):
                    l = xiser.get(j, 0)
                    if l > 0 and l <= nres1:
                        # fitzali is 0-based index, value is 1-based cd2 residue
                        fitzali[l - 1] = cd2_start + j - cd1_start

        # Convert fitzali to 0-based for wolf functions (ali[i] = j means
        # domain-local residue i aligns to cd2 residue j)
        # wolf getut/fitz expect ali[i] = j (0-based) or -1
        fitzali_0 = np.where(fitzali > 0, fitzali - 1, -1)

        # NOTE: Fortran does NOT make a fresh copy of xca between cd2 pairs
        # with the same cd1/idom. It accumulates transforms via transrotate/fitz.
        # We must replicate this behavior for ground truth matching.

        # getut: compute initial rotation/translation
        u, t, lali, rms = getut(nx, fitzali_0, xca, ca2)

        # transrotate: apply to xca (persists across cd2 pairs)
        transrotate(xca, u, t)

        # fitz: iterative refinement with initial alignment
        # (Fortran fitz takes ali as in-out, checks frozen against initial)
        fitzali_out, rms, lali, niter = fitz95(
            xca, ca2, fitzali_0, fitzrcut, fitzmaxiter)

        # Map back to full-protein alignment for scoring
        # tmpali: 1-based full residue -> 1-based cd2 residue (or 0)
        tmpali = np.zeros(nres1, dtype=np.int32)
        for i in range(nx):
            if fitzali_out[i] >= 0:
                full_res = resix[i]  # 1-based
                tmpali[full_res - 1] = fitzali_out[i] + 1  # back to 1-based

        # Score
        totscore = gettotscore95(tmpali, d1, nres1, d2, nres2)
        xscore = int(np.rint(10000 * totscore))

        # If fitz improved the score, use fitz alignment
        if xscore > score:
            # Rebuild ali from fitzali
            new_nseg = 0
            for i in range(nx):
                if fitzali_out[i] >= 0:
                    new_nseg += 1

            new_ali = [0] * (new_nseg * 4)
            j = 0
            for i in range(nx):
                if fitzali_out[i] >= 0:
                    # cd1 range: single residue (resix is 1-based)
                    new_ali[j * 2] = resix[i]
                    new_ali[j * 2 + 1] = resix[i]
                    # cd2 range: single residue (1-based)
                    new_ali[j * 2 + new_nseg * 2] = fitzali_out[i] + 1
                    new_ali[j * 2 + new_nseg * 2 + 1] = fitzali_out[i] + 1
                    j += 1

            nseg = new_nseg
            ali = new_ali
            score = xscore

        # Compute Z-score
        mean_int, sigma_int = minscore[idom]
        zscore = float(score - mean_int) / sigma_int

        # Reject low-scoring
        if zscore < zcut1:
            continue

        # Format output line (list-directed format from Fortran write(*,...))
        parts = [f'{zscore:14.7f}']
        parts.append(f'     {cd1}{cd2}')
        parts.append(f'{idom:12d}')
        parts.append(f'{score:12d}')
        parts.append(f'{nseg:12d}')
        for v in ali[:nseg*4]:
            parts.append(f'{v:12d}')
        output_lines.append(''.join(parts))

    return output_lines


def run_filter95(refine_lines, dat_dir, zcut1=1.0, fitzrcut=4.0,
                 fitzmaxiter=3):
    """
    Public API: run FILTER95 on PARSI refine output.

    Args:
        refine_lines: list of PARSI refine output line strings
        dat_dir: path to .dat files (with trailing /)
        zcut1: Z-score cutoff
        fitzrcut: FITZ distance cutoff
        fitzmaxiter: FITZ max iterations

    Returns:
        list of output line strings
    """
    return dowork_filter95(refine_lines, dat_dir, zcut1, fitzrcut, fitzmaxiter)
