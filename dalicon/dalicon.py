"""
DALICON main orchestration: Monte Carlo refinement of structural alignments.

Translated from comparemodules.f dowork_dalicon().
"""

import sys
import numpy as np
from pathlib import Path

from .scoring import gagaweights, getgagadist, RandomState, _nint
from .testi import testi
from .lean import lean_mc
from .output import produce_output

# Reuse existing dat_reader from wolf module
sys.path.insert(0, str(Path(__file__).parent.parent))
from wolf.dat_reader import read_dat
from wolf.alignment import fitz, getut
from wolf.superposition import transrotate


# Constants matching Fortran defaults
REFITOL = 2000.0
MAXITER = 100
FITZRCUT = 4.0
FITZMAXITER = 10
BLOCKSIZE = 4
IMAX0 = 50
ITRIM2 = 5
MAXRES = 5000


def _setup_new(cd, dat_dir):
    """
    Read .dat file for DALICON.

    Returns:
        nres: number of residues
        ca: (3, nres) CA coordinates (float64)
    """
    protein = read_dat(Path(dat_dir) / f'{cd}.dat')
    return protein.nres, protein.ca


def _convert_prealignment(npr, pr, nres1):
    """
    Convert prealignment ranges to per-residue array.

    The Fortran code does:
        do i=1,npr,2
            k=pr(npr*2+i)
            do j=pr(i),pr(i+1)
                preali1(j)=k; k=k+1
            end do
        end do

    Args:
        npr: number of range pairs (= nblock from input)
        pr: flat list of 4*npr integers (cd1 ranges then cd2 ranges)
        nres1: length of sequence 1

    Returns:
        preali1: (nres1+1,) int16 prealignment (1-based, 0=unaligned)
        nx: number of aligned residue pairs
    """
    preali1 = np.zeros(nres1 + 1, dtype=np.int16)
    nx = 0
    # pr is 0-based list, Fortran pr is 1-based
    # Fortran: do i=1,npr,2 → i=1,3,5,...
    # pr[i] = pr(i+1) in 0-based
    for i_fort in range(1, npr + 1, 2):
        i_py = i_fort - 1  # 0-based index into pr
        k = pr[npr * 2 + i_py]  # pr(npr*2+i) in Fortran = pr[npr*2+i-1] in 0-based
        start = pr[i_py]        # pr(i)
        end = pr[i_py + 1]      # pr(i+1)
        for j in range(start, end + 1):
            if 1 <= j <= nres1:
                preali1[j] = k
                k += 1
                nx += 1
    return preali1, nx


def dowork_dalicon(cd1, cd2, dat_dir, npr, pr, rng, lfitz=True,
                   cd1_cache=None):
    """
    Run DALICON on a single pair.

    Translated from comparemodules.f dowork_dalicon().

    Args:
        cd1, cd2: structure codes
        dat_dir: path to .dat files
        npr: number of range pairs
        pr: flat list of 4*npr ints (cd1 ranges then cd2 ranges)
        rng: RandomState for reproducible random numbers
        lfitz: whether to run fitz pre-extension
        cd1_cache: optional dict {nres, ca, d} for cd1

    Returns:
        list of WOLFITZ output lines (0 or 1)
        cd1_cache: updated cache dict
    """
    # Setup cd1
    # NOTE: In Fortran, ca1 is a module-level variable that persists across
    # comparisons with the same cd1. transrotate modifies it in-place, and
    # subsequent comparisons reuse the modified (superposed) coordinates.
    # We match this by storing the mutable ca1 array in cache and NOT copying it.
    if cd1_cache is None or cd1_cache.get('cd') != cd1:
        nres1, ca1 = _setup_new(cd1, dat_dir)
        d1 = getgagadist(ca1, nres1)
        cd1_cache = {'cd': cd1, 'nres': nres1, 'ca': ca1, 'd': d1}
    else:
        nres1 = cd1_cache['nres']
        ca1 = cd1_cache['ca']  # Same object — transrotate will modify it in-place
        d1 = cd1_cache['d']

    # Setup cd2
    nres2, ca2 = _setup_new(cd2, dat_dir)
    d2 = getgagadist(ca2, nres2)

    # Convert prealignment
    preali1, nx = _convert_prealignment(npr, pr, nres1)

    # Run fitz to extend short prealignment
    if lfitz:
        # Convert preali1 to integer array for getut (0-based)
        ali_0based = np.full(nres1, -1, dtype=np.int32)
        for i in range(1, nres1 + 1):
            if preali1[i] > 0:
                ali_0based[i - 1] = int(preali1[i]) - 1

        # Initial superposition by u3b
        u, t, nali, rms = getut(nres1, ali_0based, ca1, ca2)
        transrotate(ca1, u, t)

        # Fitz iterations
        ali_fitz, rms_fitz, nali_fitz, niter = fitz(
            ca1, ca2, FITZRCUT, FITZMAXITER)

        # Overwrite preali1 if fitzed is longer
        if nali_fitz > nx:
            for i in range(nres1):
                if ali_fitz[i] >= 0:
                    preali1[i + 1] = np.int16(ali_fitz[i] + 1)
                else:
                    preali1[i + 1] = 0

    # compresspreali is a no-op when tsil[i]=i, tsil2[j]=j

    # Testi: find tetrapeptide candidates
    ntetra, tetrapool = testi(nres1, nres2, preali1, d1, d2)

    # Initialize alignment
    ali1 = np.zeros(MAXRES + 1, dtype=np.int16)
    score = 0.0

    if ntetra > 0:
        score, ali1 = lean_mc(
            IMAX0, BLOCKSIZE, nres1, nres2, d1, d2, True,
            preali1, ITRIM2, ntetra, tetrapool, rng)

    # Iterate until improvement < refitol
    oldscore = score
    ds = oldscore
    iter_count = 0

    while ds > REFITOL and iter_count < MAXITER:
        iter_count += 1
        # Copy ali1 to preali1
        for i in range(1, nres1 + 1):
            preali1[i] = ali1[i]

        ntetra, tetrapool = testi(nres1, nres2, preali1, d1, d2)
        ali1_new = np.zeros(MAXRES + 1, dtype=np.int16)
        score_new = 0.0

        if ntetra > 0:
            score_new, ali1_new = lean_mc(
                IMAX0, BLOCKSIZE, nres1, nres2, d1, d2, True,
                preali1, ITRIM2, ntetra, tetrapool, rng)

        ali1 = ali1_new
        score = score_new
        ds = score - oldscore
        oldscore = score

    # Produce output
    lines = produce_output(cd1, cd2, nres1, nres2, ali1, score)

    return lines, cd1_cache


def parse_dalicon_input(filepath):
    """
    Parse DALICON input file (grouped by cd1 with END delimiters).

    Returns:
        records: list of (cd1, cd2, nblock, data_values) tuples
    """
    records = []
    with open(filepath) as f:
        lines = [l.rstrip() for l in f.readlines()]

    i = 0
    cd1 = None
    while i < len(lines):
        line = lines[i].strip()

        # END marks boundary between cd1 groups
        if line == 'END':
            i += 1
            # Next non-empty line is new cd1
            while i < len(lines) and not lines[i].strip():
                i += 1
            if i < len(lines) and lines[i].strip() != 'END':
                cd1 = lines[i].strip()
                i += 1
            continue

        if not line:
            i += 1
            continue

        if cd1 is None:
            cd1 = line
            i += 1
            continue

        # cd2 name (possibly with * suffix)
        if not line.replace('-', '').replace('.', '').isdigit():
            cd2 = line.rstrip('*')
            i += 1
            # Skip blank lines
            while i < len(lines) and not lines[i].strip():
                i += 1
            if i >= len(lines):
                break
            nblock = int(lines[i].strip())
            i += 1

            # Read 4*nblock integer values
            values = []
            while len(values) < nblock * 4 and i < len(lines):
                vline = lines[i].strip()
                if not vline or vline == 'END':
                    break
                try:
                    values.extend([int(x) for x in vline.split()])
                except ValueError:
                    break
                i += 1

            records.append((cd1, cd2, nblock, values))
        else:
            i += 1

    return records


def run_dalicon(input_file, dat_dir):
    """
    Run DALICON on all pairs from input file.

    Args:
        input_file: path to DALICON input file
        dat_dir: directory containing .dat files

    Returns:
        output_lines: list of WOLFITZ output strings
    """
    # Initialize weights
    gagaweights()

    # Create random state (seed persists across all comparisons)
    rng = RandomState(seed=124462287)

    records = parse_dalicon_input(input_file)

    output_lines = []
    cd1_cache = None

    for cd1, cd2, nblock, values in records:
        lines, cd1_cache = dowork_dalicon(
            cd1, cd2, dat_dir, nblock, values, rng,
            lfitz=True, cd1_cache=cd1_cache)
        output_lines.extend(lines)

    return output_lines
