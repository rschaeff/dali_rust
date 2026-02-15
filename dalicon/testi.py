"""
Tetrapeptide seeding for DALICON Monte Carlo.

Translated from testi.f and ssap.f (getleftrite).
"""

import numpy as np
from .scoring import scorefun

MAXPAIR = 6000


def getleftrite(nres1, nres2, preali1, width=10):
    """
    Compute left/right bounds for search space based on prealignment.

    Translated from ssap.f getleftrite(). Uses 1-based indexing.

    Args:
        nres1: length of sequence 1
        nres2: length of sequence 2
        preali1: (nres1+1,) prealignment array (1-based, 0=unaligned)
        width: half-width of search band

    Returns:
        left: (nres1+1,) left bounds (1-based)
        rite: (nres1+1,) right bounds (1-based)
    """
    left = np.zeros(nres1 + 1, dtype=np.int32)
    rite = np.full(nres1 + 1, nres2, dtype=np.int32)

    # Forward pass
    for i in range(2, nres1 + 1):
        if preali1[i] > 0:
            left[i] = max(left[i - 1], preali1[i] - width)
        else:
            left[i] = left[i - 1]

    # Backward pass
    for i in range(nres1 - 1, 0, -1):
        if preali1[i] > 0:
            rite[i] = min(preali1[i] + width, rite[i + 1])
        else:
            rite[i] = rite[i + 1]

    return left, rite


def testi(nres1, nres2, preali1, d1, d2):
    """
    Build tetrapeptide candidate pool from prealignment and distance maps.

    Translated from testi.f testi().

    Args:
        nres1, nres2: sequence lengths
        preali1: (nres1+1,) prealignment (1-based, 0=unaligned)
        d1: (nres1+1, nres1+1) int16 distance matrix (1-based)
        d2: (nres2+1, nres2+1) int16 distance matrix (1-based)

    Returns:
        ntetra: number of tetrapeptide candidates
        tetrapool: list of (i, j) tuples (1-based residue indices)
    """
    # Collect aligned positions
    aligned = []
    for i in range(1, nres1 + 1):
        if preali1[i] != 0:
            aligned.append(i)

    left, rite = getleftrite(nres1, nres2, preali1, 10)

    # Build score map
    # Use a dictionary for sparse storage (map can be huge)
    # map[i,j] = sum of scorefun for all aligned pairs
    score_map = {}
    for i in range(1, nres1 + 1):
        for j in range(left[i] + 1, rite[i]):
            if j < 1 or j > nres2:
                continue
            x = 0.0
            for i0 in aligned:
                j0 = preali1[i0]
                if i != i0 and j != j0:
                    x += scorefun(d1[i0, i], d2[j0, j])
            score_map[(i, j)] = x

    # Accept all bands +-3 residues around aligned segments
    for i in range(1, nres1 + 1):
        k = preali1[i]
        if k > 0:
            for j in range(max(1, k - 3), min(nres1 - 3, k + 3) + 1):
                score_map[(i, j)] = 1.0

    # Count tetrapeptides (4 consecutive positive map values)
    tetrapool = []
    for i in range(1, nres1 - 2):  # i goes 1..nres1-3
        for j in range(left[i] + 1, rite[i]):
            if j < 1 or j > nres2 - 3:
                continue
            x = 0.0
            for l in range(4):
                x += score_map.get((i + l, j + l), 0.0)
            if x > 0.0 and len(tetrapool) < MAXPAIR:
                tetrapool.append((i, j))

    return len(tetrapool), tetrapool
