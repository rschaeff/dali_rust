"""
DP module: Z-score computation from WOLFITZ alignments.

Translated from comparemodules.f dp module.
"""

import numpy as np
from dataclasses import dataclass
from pathlib import Path

# Constants
D0 = 0.20
ENVELOPE_RADIUS = 20.0
ZCUT_DEFAULT = 2.0


def _nint(x):
    """Fortran-compatible rounding: round half away from zero.

    Fortran's nint(0.5)=1, nint(1.5)=2, nint(2.5)=3.
    Python's round/np.rint use banker's rounding (half to even).
    """
    if x >= 0:
        return int(x + 0.5)
    else:
        return int(x - 0.5)


@dataclass
class DCCPResult:
    cd1: str
    cd2: str
    score: float   # x1: score from idom1=1, idom2=1
    zmax: float    # max z-score across all domain pairs
    rmsd: float    # placeholder (9.9)
    lali: int      # alignment length (0 — getide commented out in Fortran)
    ide: int       # sequence identity (0 — getide commented out)
    nblock: int
    l1: list       # start positions (1-based)
    r1: list       # end positions (1-based)
    l2: list       # start positions (1-based)
    r2: list       # end positions (1-based)

    def format_dccp(self):
        """Format as DCCP output matching Fortran format 601/605/610."""
        lines = []
        # Format 601: (1x,a9,f8.1,f4.1,i4,f8.1,8x,i4,3x,i4,16x,a5,1x,a5)
        header = (f" DCCP   1 {self.score:8.1f}{self.rmsd:4.1f}{self.lali:4d}"
                  f"{self.zmax:8.1f}{'':8s}{self.ide:4d}{'':3s}{self.nblock:4d}"
                  f"{'':16s}{self.cd1:>5s} {self.cd2:>5s}")
        lines.append(header)
        # Format 605: (1x,a9)
        lines.append(' alignment')
        # Format 610: (8(i4,2x,i4)) — cd1 ranges
        lines.extend(_format_ranges(self.l1, self.r1, self.nblock))
        # Format 610: (8(i4,2x,i4)) — cd2 ranges
        lines.extend(_format_ranges(self.l2, self.r2, self.nblock))
        return '\n'.join(lines)


def _format_ranges(l, r, nblock):
    """Format alignment ranges as Fortran (8(i4,2x,i4)) output."""
    range_lines = []
    vals = []
    for b in range(nblock):
        vals.append(f"{l[b]:4d}{'':2s}{r[b]:4d}")
    # 8 groups per line
    for i in range(0, len(vals), 8):
        range_lines.append(''.join(vals[i:i + 8]))
    return range_lines


def dpweights():
    """Compute Gaussian envelope weights.

    wght[i] = exp(-(i/radius)^2) for i=0..100

    Returns:
        wght: numpy array of shape (101,)
    """
    wght = np.zeros(101, dtype=np.float64)
    x = 1.0 / (ENVELOPE_RADIUS * ENVELOPE_RADIUS)
    for i in range(101):
        wght[i] = np.exp(-x * i * i)
    return wght


def dpgetdist(ca, nres):
    """Compute pairwise CA distance matrix as integer*2 (dist * 10).

    Translated from comparemodules.f dpgetdist().

    Args:
        ca: (3, nres) CA coordinates
        nres: number of residues

    Returns:
        d: (nres, nres) int16 distance matrix (distances * 10, rounded)
    """
    # Use float32 to match Fortran's single-precision arithmetic
    ca32 = ca.astype(np.float32) if ca.dtype != np.float32 else ca
    d = np.zeros((nres, nres), dtype=np.int16)
    ten = np.float32(10.0)
    for i in range(nres):
        for j in range(i):
            dx = ca32[0, i] - ca32[0, j]
            dy = ca32[1, i] - ca32[1, j]
            dz = ca32[2, i] - ca32[2, j]
            dist = np.sqrt(dx * dx + dy * dy + dz * dz)
            x = _nint(float(ten * dist))
            # Clip to int16 range
            if x > 32767:
                x = 32767
            d[i, j] = x
            d[j, i] = x
    return d


def dpscorefun(a, b, wght):
    """Distance-based scoring function.

    Translated from comparemodules.f dpscorefun().

    Args:
        a, b: int16 distances (dist * 10)
        wght: weight array from dpweights()

    Returns:
        s: score value
    """
    x = abs(int(a) - int(b)) / 10.0
    y = (int(a) + int(b)) / 20.0
    if y > 100.0:
        return 0.0
    iy = _nint(y)
    if iy < 0:
        iy = 0
    if iy > 100:
        iy = 100
    if y > 0.0:
        return float(wght[iy]) * (D0 - x / y)
    else:
        return float(wght[iy]) * D0


def totscore(ali1, nres1, d1, d2, wght):
    """Compute total alignment score.

    Double loop over all pairs of aligned residues.

    Translated from comparemodules.f totscore().

    Args:
        ali1: (nres1,) alignment array. ali1[i]=j (1-based) means
              residue i+1 aligned to residue j. 0 means unaligned.
        nres1: number of residues in protein 1
        d1: (nres1, nres1) distance matrix for protein 1
        d2: (nres2, nres2) distance matrix for protein 2
        wght: weight array

    Returns:
        tots: total score
    """
    # Collect aligned residue indices (0-based)
    a = []
    for i in range(nres1):
        if ali1[i] != 0:
            a.append(i)

    tots = 0.0
    n = len(a)
    for i in range(n):
        k = a[i]
        for j in range(n):
            l = a[j]
            q = abs(ali1[k]) - 1  # Convert 1-based target to 0-based
            r = abs(ali1[l]) - 1
            tots += dpscorefun(d1[k, l], d2[q, r], wght)

    return tots


def zscore_func(l1, l2, score):
    """Compute Z-score from domain lengths and raw score.

    Uses cubic polynomial calibration.

    Translated from comparemodules.f zscore().

    Args:
        l1: length of domain 1 (number of residues)
        l2: length of domain 2
        score: raw alignment score

    Returns:
        z: Z-score
    """
    n12 = np.sqrt(float(l1 * l2))
    x = min(n12, 400.0)
    mean = 7.9494 + 0.70852 * x + 2.5895e-4 * x * x - 1.9156e-6 * x * x * x
    if n12 > 400.0:
        mean = mean + (n12 - 400.0) * 1.0
    sigma = 0.50 * mean
    z = (score - mean) / max(1.0, sigma)
    return z


def dpsetup(filepath):
    """Read .dat file for DP: extract CA coords and domain definitions.

    Reads three >>>> sections:
    - 1st: nres, nseg, skip segment ranges, read CA coords
    - 2nd: skip (SSE hierarchy)
    - 3rd: domain definitions

    Filters domains to keep only first domain, '+', or '*' types.

    Translated from comparemodules.f dpsetup().

    Args:
        filepath: path to .dat file

    Returns:
        nres: number of residues
        ca: (3, nres) CA coordinates
        d: (nres, nres) int16 distance matrix
        ndom: number of domains (after filtering)
        domns: list of segment counts per domain
        domseglist: list of list of (start, end) per domain (1-based)
    """
    with open(filepath, errors='replace') as f:
        lines = f.readlines()

    idx = 0
    header_count = 0
    nres = 0
    ca = None

    # Domain info from 3rd header
    domns_raw = {}
    domseglist_raw = {}
    node_types = {}
    max_j = 0

    while idx < len(lines):
        line = lines[idx]
        if not line.startswith('>>>>'):
            idx += 1
            continue

        header_count += 1

        if header_count == 1:
            # First >>>>: read nres and nseg, then CA coords
            # Format: (10x,i5) for ndom, but we read nres/nseg from chars 11-20
            fields = line[10:].split()
            nres = int(fields[0])
            nseg = int(fields[1])
            idx += 1

            # Skip nseg segment range lines
            for _ in range(nseg):
                idx += 1

            # Read CA coordinates: free-format, 3*nres values
            ca_values = []
            total_needed = 3 * nres
            while len(ca_values) < total_needed and idx < len(lines):
                vals = lines[idx].split()
                ca_values.extend(float(v) for v in vals)
                idx += 1

            ca = np.zeros((3, nres), dtype=np.float64)
            for j in range(nres):
                ca[0, j] = ca_values[j * 3]
                ca[1, j] = ca_values[j * 3 + 1]
                ca[2, j] = ca_values[j * 3 + 2]
            continue

        elif header_count == 2:
            # Second >>>>: skip
            idx += 1
            continue

        elif header_count == 3:
            # Third >>>>: read domain definitions
            # ndom from format (10x,i5)
            ndom_header = int(line[10:15])
            idx += 1

            # Read ndom_header domain lines with format 710:
            # (i4,1x,a1,13x,400i4)
            for _ in range(ndom_header):
                if idx >= len(lines):
                    break
                dline = lines[idx]
                idx += 1

                if len(dline.rstrip()) < 6:
                    continue

                # Parse: i4 (pos 0-3), 1x (pos 4), a1 (pos 5), 13x (pos 6-18)
                j = int(dline[0:4])
                node = dline[5]

                # 400i4 starting at position 19
                int_str = dline[19:]
                ints = []
                pos = 0
                while pos + 4 <= len(int_str.rstrip('\n')):
                    chunk = int_str[pos:pos + 4].strip()
                    if chunk:
                        try:
                            ints.append(int(chunk))
                        except ValueError:
                            break
                    pos += 4

                if ints:
                    ns = ints[0]  # domns(j)
                    segs = []
                    for s in range(ns):
                        if 1 + s * 2 + 1 < len(ints):
                            start = ints[1 + s * 2]
                            end = ints[1 + s * 2 + 1]
                            segs.append((start, end))
                    domns_raw[j] = ns
                    domseglist_raw[j] = segs
                else:
                    domns_raw[j] = 0
                    domseglist_raw[j] = []

                node_types[j] = node
                max_j = j
            break

        idx += 1

    # Compute distance matrix
    d = dpgetdist(ca, nres)

    # Filter domains: keep i==1 or node_type in ('+', '*')
    # Matches Fortran: if(i.eq.1.or.node_type(i).eq.'+'.or.node_type(i).eq.'*')
    ndom = 0
    domns = []
    domseglist = []
    for i in range(1, max_j + 1):
        if i not in node_types:
            continue
        if i == 1 or node_types[i] == '+' or node_types[i] == '*':
            domns.append(domns_raw.get(i, 0))
            domseglist.append(domseglist_raw.get(i, []))
            ndom += 1

    return nres, ca, d, ndom, domns, domseglist


def dopair(idom1, idom2, l1, r1, l2, r2, nblock, nres1, nres2,
           d1, d2, domns1, domseglist1, domns2, domseglist2, wght):
    """Score a domain pair.

    Translated from comparemodules.f dopair().

    Args:
        idom1, idom2: domain indices (0-based)
        l1, r1, l2, r2: alignment block ranges (1-based, lists)
        nblock: number of blocks
        nres1, nres2: residue counts
        d1, d2: distance matrices
        domns1, domseglist1: domain definitions for protein 1
        domns2, domseglist2: domain definitions for protein 2
        wght: weight array

    Returns:
        z: Z-score for this domain pair
        x: raw score for this domain pair
    """
    # Mark active residues in each domain (1-based indexing)
    lactive1 = [False] * (nres1 + 1)
    lactive2 = [False] * (nres2 + 1)

    len1 = 0
    for seg_start, seg_end in domseglist1[idom1]:
        for j in range(seg_start, min(seg_end + 1, nres1 + 1)):
            lactive1[j] = True
            len1 += 1

    len2 = 0
    for seg_start, seg_end in domseglist2[idom2]:
        for j in range(seg_start, min(seg_end + 1, nres2 + 1)):
            lactive2[j] = True
            len2 += 1

    # Construct ali1 (1-based: ali1[i] = j means residue i aligned to j)
    ali1 = np.zeros(nres1 + 1, dtype=np.int32)
    lali = 0
    for b in range(nblock):
        for j in range(l1[b], r1[b] + 1):
            j2 = l2[b] + j - l1[b]
            if 1 <= j <= nres1 and 1 <= j2 <= nres2:
                if lactive1[j] and lactive2[j2]:
                    ali1[j] = j2
                    lali += 1

    if lali == 0:
        return 0.0, 0.0

    # Convert to 0-based for totscore
    # ali1_0[i] = ali1[i+1] (keeping 1-based target values)
    ali1_0 = np.zeros(nres1, dtype=np.int32)
    for i in range(1, nres1 + 1):
        ali1_0[i - 1] = ali1[i]

    x = totscore(ali1_0, nres1, d1, d2, wght)
    z = zscore_func(len1, len2, x)

    return z, x


def compressblocks(nblock, l1, r1, l2, r2):
    """Merge consecutive single-residue alignment blocks.

    Exact translation of comparemodules.f compressblocks().
    """
    if nblock == 0:
        return 0

    i = 0
    n = 0

    while i < nblock:
        if n < len(l1):
            l1[n] = l1[i]
            l2[n] = l2[i]
            r1[n] = r1[i]
            r2[n] = r2[i]
        else:
            l1.append(l1[i])
            l2.append(l2[i])
            r1.append(r1[i])
            r2.append(r2[i])

        if i < nblock - 1:
            while (l1[i] == r1[i] and l2[i] == r2[i] and
                   l1[i + 1] == r1[i] + 1 and l2[i + 1] == r2[i] + 1):
                if i + 1 >= nblock - 1:
                    break
                i += 1
            r1[n] = r1[i]
            r2[n] = r2[i]

        n += 1
        i += 1

    return n


def run_dp(wolfitz_input, dat_dir, zcut=ZCUT_DEFAULT):
    """Run DP scoring on WOLFITZ input lines.

    Mimics the Fortran dowork_dp loop: for each WOLFITZ hit, setup proteins,
    score all domain pairs, output DCCP if zmax >= zcut.

    Args:
        wolfitz_input: path to WOLFITZ output file
        dat_dir: directory containing .dat files
        zcut: Z-score cutoff (default 2.0)

    Returns:
        list of DCCPResult
    """
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from validation.parsers import parse_wolfitz

    dat_dir = Path(dat_dir)
    hits = parse_wolfitz(wolfitz_input)

    # Compute weights once
    wght = dpweights()

    # Cache protein data
    protein_cache = {}

    results = []

    for hit in hits:
        cd1 = hit.cd1
        cd2 = hit.cd2
        nblock = hit.nblock

        # Parse ranges from raw_values
        l1 = []
        r1 = []
        l2 = []
        r2 = []
        for b in range(nblock):
            l1.append(hit.raw_values[b * 2])
            r1.append(hit.raw_values[b * 2 + 1])
        offset = nblock * 2
        for b in range(nblock):
            l2.append(hit.raw_values[offset + b * 2])
            r2.append(hit.raw_values[offset + b * 2 + 1])

        # Setup proteins (cached)
        if cd1 not in protein_cache:
            protein_cache[cd1] = dpsetup(dat_dir / f'{cd1}.dat')
        if cd2 not in protein_cache:
            protein_cache[cd2] = dpsetup(dat_dir / f'{cd2}.dat')

        nres1, ca1, d1, ndom1, domns1, domseglist1 = protein_cache[cd1]
        nres2, ca2, d2, ndom2, domns2, domseglist2 = protein_cache[cd2]

        if ndom1 == 0 or ndom2 == 0:
            continue

        # Score all domain pairs
        zmax = 0.0
        x1 = 0.0
        rmsd = 9.9

        for idom1 in range(ndom1):
            for idom2 in range(ndom2):
                z, x = dopair(idom1, idom2, l1, r1, l2, r2, nblock,
                              nres1, nres2, d1, d2,
                              domns1, domseglist1, domns2, domseglist2,
                              wght)
                if z > zmax:
                    zmax = z
                if idom1 == 0 and idom2 == 0:
                    x1 = x

        # Output if any domain pair passed z cutoff
        if zmax >= zcut:
            # Apply compressblocks (should be idempotent on WOLF output)
            l1_c = l1[:]
            r1_c = r1[:]
            l2_c = l2[:]
            r2_c = r2[:]
            nblock_c = compressblocks(nblock, l1_c, r1_c, l2_c, r2_c)
            l1_c = l1_c[:nblock_c]
            r1_c = r1_c[:nblock_c]
            l2_c = l2_c[:nblock_c]
            r2_c = r2_c[:nblock_c]

            results.append(DCCPResult(
                cd1=cd1, cd2=cd2,
                score=x1, zmax=zmax,
                rmsd=rmsd, lali=0, ide=0,
                nblock=nblock_c,
                l1=l1_c, r1=r1_c, l2=l2_c, r2=r2_c,
            ))

    return results
