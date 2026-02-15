"""
WOLF top-level orchestration: compare two proteins via SSE hashing.

Translated from comparemodules.f dowork_wolf() and wolf_original.f.
"""

import numpy as np
from dataclasses import dataclass
from pathlib import Path

from .dat_reader import read_dat, Protein
from .geometry import vec, twist, preparex, compute_neidist
from .spatial_hash import SpatialHashGrid, fung, lgrid
from .alignment import fitz


# Constants matching Fortran serialcompare.f
RCUT_FITZ = 4.0
MAXITER = 20
NEIBORCUTOFF = 12.0


@dataclass
class WolfResult:
    cd1: str
    cd2: str
    nblock: int
    l1: list  # start positions in cd1 (1-based)
    r1: list  # end positions in cd1 (1-based)
    l2: list  # start positions in cd2 (1-based)
    r2: list  # end positions in cd2 (1-based)

    def to_wolfitz_line(self):
        """Format as WOLFITZ output line matching Fortran free-format write."""
        # Fortran: write(*,*) 'WOLFITZ ',cd1,cd2,nblock,(l1(i),r1(i),...),(l2(i),r2(i),...)
        # cd1 and cd2 are character*5, space-padded
        cd1_padded = f'{self.cd1:<5s}'
        cd2_padded = f'{self.cd2:<5s}'
        parts = [f'WOLFITZ {cd1_padded}{cd2_padded}']
        # nblock and all values as right-justified integers (Fortran list-directed)
        vals = [self.nblock]
        for i in range(self.nblock):
            vals.append(self.l1[i])
            vals.append(self.r1[i])
        for i in range(self.nblock):
            vals.append(self.l2[i])
            vals.append(self.r2[i])
        parts.extend(f'{v:>12d}' for v in vals)
        return ''.join(parts)


def _distance(a, b):
    """Euclidean distance between two 3D points."""
    d = a - b
    return np.sqrt(np.sum(d * d))


def _cosi(z, x):
    """Cosine of angle between two 3D vectors."""
    lz = np.sqrt(np.sum(z * z))
    lx = np.sqrt(np.sum(x * x))
    norm = max(1e-6, lz * lx)
    return np.sum(z * x) / norm


def setup_protein(filepath):
    """
    Read .dat file and compute SSE vectors and neighbor distances.

    Returns:
        protein: Protein dataclass
        midpoint: (3, nseg)
        direction: (3, nseg)
        neidist: (nseg, nseg)
    """
    protein = read_dat(filepath)
    midpoint, direction = vec(protein.ca, protein.nres, protein.nseg,
                              protein.segmentrange)
    neidist = compute_neidist(protein.nseg, midpoint)
    return protein, midpoint, direction, neidist


def load_protein(grid, nseg, midpoint, direction, neidist, rcut):
    """
    Build spatial hash grid for a protein's SSE descriptors.

    For each SSE pair (iseg, jseg) within distance cutoff:
      - Build canonical frame via preparex + twist
      - Insert all other SSEs into grid via boxit

    Args:
        grid: SpatialHashGrid to populate
        nseg: number of SSEs
        midpoint: (3, nseg)
        direction: (3, nseg)
        neidist: (nseg, nseg) pairwise distances
        rcut: distance cutoff for SSE pairs (neiborcutoff)
    """
    grid.clear()
    ns = 0
    for iseg in range(nseg):
        for jseg in range(nseg):
            if jseg != iseg and neidist[iseg, jseg] < rcut:
                x = preparex(iseg, jseg, nseg, midpoint, direction)
                twist(x, 3 + 2 * nseg)
                grid.boxit(nseg, x, iseg, jseg)
                ns += 1
    return ns


def compare(nseg_cd2, mid_cd2, dir_cd2, secstr_cd2,
            nseg_cd1, secstr_cd1,
            grid, rcut, neidist_cd2):
    """
    Find best SSE pair mapping between cd2 (query) and cd1 (in grid).

    For each cd2 SSE pair, transform cd2's other SSEs into canonical frame,
    look up matching cd1 entries in the spatial hash grid, accumulate votes.

    Args:
        nseg_cd2: cd2's SSE count
        mid_cd2, dir_cd2: cd2's SSE midpoints and directions
        secstr_cd2: cd2's SSE types
        nseg_cd1: used as protein_nseg (for vote counting bounds)
        secstr_cd1: cd1's SSE types
        grid: SpatialHashGrid containing cd1's SSE data
        rcut: distance cutoff (neiborcutoff)
        neidist_cd2: cd2's pairwise SSE distances

    Returns:
        bestpair: (iseg, jseg, useg, vseg) — cd2 pair (iseg, jseg), cd1 pair (useg, vseg)
                  All 0-based. Returns None if no match found.
        protcount: best vote count
    """
    protcount = 0
    bestpair = None

    for iseg in range(nseg_cd2):
        atype = secstr_cd2[iseg]
        for jseg in range(nseg_cd2):
            lest = (iseg < jseg)  # topology flag for iseg vs jseg ordering
            if iseg == jseg:
                continue
            r = neidist_cd2[iseg, jseg]
            if r >= rcut:
                continue

            # Initialize vote counts
            # Fortran uses count(maxseg, maxseg), initialized to 0 up to protein_nseg
            # We use a dict for efficiency
            count = {}

            # Build canonical frame for cd2 pair (iseg, jseg)
            x = preparex(iseg, jseg, nseg_cd2, mid_cd2, dir_cd2)
            twist(x, 3 + 2 * nseg_cd2)

            # For each other cd2 SSE, look up in cd1's grid
            for kseg in range(nseg_cd2):
                if kseg == iseg:
                    continue
                less = (iseg < kseg)  # topology flag for iseg vs kseg ordering

                # Compute midpoint and direction of kseg in canonical frame
                midx = (x[:, 3 + kseg] + x[:, 3 + nseg_cd2 + kseg]) / 2.0
                dirx = x[:, 3 + nseg_cd2 + kseg] - midx

                ctype = secstr_cd2[kseg]

                # Grid position
                gx0 = fung(midx[0])
                gy0 = fung(midx[1])
                gz0 = fung(midx[2])

                # Search grid neighborhood
                for entry in grid.lookup(gx0, gy0, gz0, radius=2):
                    cur_a, cur_b, cur_c, link_from, link_to = entry

                    # Filter: secondary structure type match
                    if secstr_cd1[cur_a] != atype:
                        continue
                    if secstr_cd1[cur_c] != ctype:
                        continue

                    # Topology filters
                    if less:
                        if cur_a >= cur_c:
                            continue
                    else:
                        if cur_a <= cur_c:
                            continue

                    if lest:
                        if cur_a >= cur_b:
                            continue
                    else:
                        if cur_a <= cur_b:
                            continue

                    # Midpoint distance filter (< 4.0 Å)
                    mid_entry = (link_from + link_to) / 2.0
                    dir_entry = link_to - mid_entry
                    if _distance(mid_entry, midx) > 4.0:
                        continue

                    # Direction cosine filter (> 0.5)
                    if _cosi(dir_entry, dirx) < 0.5:
                        continue

                    # Vote
                    key = (cur_a, cur_b)
                    count[key] = count.get(key, 0) + 1

            # Find best (useg, vseg) for this cd2 pair — two-phase search
            # matching Fortran exactly:
            # Phase 1: find local best within protein_nseg bounds
            # Phase 2: compare to global best
            # The Fortran iterates lseg=1..protein_nseg, kseg=1..protein_nseg
            # (1-based), so we iterate 0..protein_nseg-1 (0-based).
            protein_nseg = nseg_cd1  # cd2's nseg passed as protein_nseg
            local_best = 0
            local_useg = -1
            local_vseg = -1
            for lseg in range(protein_nseg):
                for kseg in range(protein_nseg):
                    c = count.get((lseg, kseg), 0)
                    if c > local_best:
                        local_best = c
                        local_useg = lseg
                        local_vseg = kseg

            if local_best > protcount:
                protcount = local_best
                bestpair = (iseg, jseg, local_useg, local_vseg)

    return bestpair, protcount


def compressblocks(nblock, l1, r1, l2, r2):
    """
    Merge consecutive single-residue alignment blocks.

    Exact translation of comparemodules.f compressblocks().

    Args:
        nblock: number of blocks
        l1, r1, l2, r2: lists of block boundaries (1-based, modified in-place)

    Returns:
        new_nblock: compressed block count
    """
    if nblock == 0:
        return 0

    i = 0  # read pointer (0-based)
    n = 0  # write pointer (0-based)

    while i < nblock:
        # Copy current block
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
            # Try to merge consecutive single-residue blocks
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


def wolf_compare(cd1, cd2, dat_dir, grid=None, cd1_data=None):
    """
    Compare two proteins using the WOLF SSE hashing algorithm.

    Args:
        cd1: query protein code (e.g. "101mA")
        cd2: target protein code (e.g. "1a00A")
        dat_dir: directory containing .dat files
        grid: optional pre-built SpatialHashGrid for cd1
        cd1_data: optional pre-computed cd1 data tuple
                  (protein, midpoint, direction, neidist)

    Returns:
        WolfResult or None if no alignment found.
    """
    dat_dir = Path(dat_dir)

    # Setup cd1
    if cd1_data is None:
        prot1, mid1, dir1, neidist1 = setup_protein(dat_dir / f'{cd1}.dat')
    else:
        prot1, mid1, dir1, neidist1 = cd1_data

    # Check cd1 has enough SSEs
    if prot1.nseg <= 2:
        return None

    # Build grid for cd1 if not provided
    if grid is None:
        grid = SpatialHashGrid()
        load_protein(grid, prot1.nseg, mid1, dir1, neidist1, NEIBORCUTOFF)

    # Setup cd2
    prot2, mid2, dir2, neidist2 = setup_protein(dat_dir / f'{cd2}.dat')

    # Check cd2 has enough SSEs
    if prot2.nseg <= 2:
        return None

    # Compare: cd2 SSE pairs searched against cd1's grid
    # Note: Fortran passes nseg (cd2's) as protein_nseg
    bp, protcount = compare(
        prot2.nseg, mid2, dir2, prot2.secstr,
        prot2.nseg, prot1.secstr,  # protein_nseg = cd2's nseg (Fortran behavior)
        grid, NEIBORCUTOFF, neidist2
    )

    if bp is None or bp[0] == bp[1]:
        return None
    if bp[2] == bp[3]:
        return None

    iseg, jseg, useg, vseg = bp

    # Build initial superposition for cd1 (using bestpair SSEs from cd1)
    x = preparex(useg, vseg, 0, mid1, dir1)
    # Append cd1's CA coords after the 3 reference points
    x_full = np.zeros((3, 3 + prot1.nres), dtype=np.float64)
    x_full[:, :3] = x[:, :3]
    for j in range(prot1.nres):
        x_full[:, j + 3] = prot1.ca[:, j]
    twist(x_full, 3 + prot1.nres)
    # Extract transformed CA coords (skip 3 reference points)
    xca = x_full[:, 3:3 + prot1.nres].copy()

    # Build initial superposition for cd2 (using bestpair SSEs from cd2)
    y = preparex(iseg, jseg, 0, mid2, dir2)
    y_full = np.zeros((3, 3 + prot2.nres), dtype=np.float64)
    y_full[:, :3] = y[:, :3]
    for j in range(prot2.nres):
        y_full[:, j + 3] = prot2.ca[:, j]
    twist(y_full, 3 + prot2.nres)
    yca = y_full[:, 3:3 + prot2.nres].copy()

    # Iterative alignment refinement
    ali, rms, lali, niter = fitz(xca, yca, RCUT_FITZ, MAXITER)

    # Convert alignment to blocks (1-based residue indices)
    l1, r1, l2, r2 = [], [], [], []
    nblock = 0
    for i in range(prot1.nres):
        if ali[i] >= 0:
            nblock += 1
            l1.append(i + 1)   # convert to 1-based
            r1.append(i + 1)
            l2.append(ali[i] + 1)
            r2.append(ali[i] + 1)

    if nblock == 0:
        return None

    # Compress blocks
    nblock = compressblocks(nblock, l1, r1, l2, r2)
    l1 = l1[:nblock]
    r1 = r1[:nblock]
    l2 = l2[:nblock]
    r2 = r2[:nblock]

    return WolfResult(
        cd1=cd1, cd2=cd2, nblock=nblock,
        l1=l1, r1=r1, l2=l2, r2=r2,
    )


def run_wolf_all_pairs(structures, dat_dir):
    """
    Run WOLF comparison for all pairs of structures.

    Mimics the Fortran serialcompare loop: for each cd1, build grid once,
    then compare against all cd2s.

    Args:
        structures: list of structure codes
        dat_dir: directory containing .dat files

    Returns:
        results: list of WolfResult (only non-None results)
    """
    dat_dir = Path(dat_dir)
    results = []

    for cd1 in structures:
        # Setup cd1 once
        cd1_data = setup_protein(dat_dir / f'{cd1}.dat')
        prot1, mid1, dir1, neidist1 = cd1_data

        if prot1.nseg <= 2:
            continue

        # Build grid for cd1
        grid = SpatialHashGrid()
        load_protein(grid, prot1.nseg, mid1, dir1, neidist1, NEIBORCUTOFF)

        for cd2 in structures:
            result = wolf_compare(cd1, cd2, dat_dir, grid=grid, cd1_data=cd1_data)
            if result is not None:
                results.append(result)

    return results
