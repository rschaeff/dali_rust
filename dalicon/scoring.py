"""
Scoring functions and distance matrices for DALICON.

Translated from gagatool.f: scorefun, gagaweights, getgagadist, gagadistance.
Also ran() from ran.f and getp() from gagatool.f.
"""

import numpy as np
import math


# Weight table (Gaussian envelope), filled by gagaweights()
# Use float32 to match Fortran real precision
_weight = np.zeros(101, dtype=np.float32)


def _nint(x):
    """Fortran-compatible nearest integer (rounds 0.5 away from zero)."""
    if x >= 0:
        return int(x + 0.5)
    else:
        return -int(-x + 0.5)


def gagaweights():
    """Compute Gaussian weight table: weight[i] = exp(-0.0025 * i^2)."""
    global _weight
    for i in range(101):
        # Use float32 to match Fortran real precision
        _weight[i] = np.float32(math.exp(np.float32(-0.0025) * np.float32(i * i)))


def scorefun(a, b):
    """
    Elastic distance scoring function.

    Translated from gagatool.f scorefun().
    Uses float32 arithmetic to match Fortran real precision.

    Args:
        a, b: int16 distance values (in 0.1 Angstrom units)

    Returns:
        score: float32
    """
    ia = int(a)
    ib = int(b)
    x = np.float32(abs(ia - ib)) / np.float32(10.0)
    y = np.float32(ia + ib) / np.float32(20.0)
    yint = _nint(float(y))
    if y > np.float32(100.0):
        return np.float32(0.0)
    if y > np.float32(0.0):
        s = np.float32(_weight[yint]) * (np.float32(0.20) - x / y)
    else:
        s = np.float32(_weight[yint]) * np.float32(0.20)
    return np.float32(s * np.float32(100.0))


def gagadistance(ca, i, j):
    """
    Compute integer distance between two CA atoms (in 0.1 A units).

    Args:
        ca: (3, nres) coordinate array
        i, j: residue indices (0-based)

    Returns:
        int16 distance
    """
    dx = ca[0, i] - ca[0, j]
    dy = ca[1, i] - ca[1, j]
    dz = ca[2, i] - ca[2, j]
    return np.int16(_nint(10.0 * math.sqrt(dx*dx + dy*dy + dz*dz)))


def getgagadist(ca, nres):
    """
    Compute full integer distance matrix.

    Args:
        ca: (3, nres) CA coordinates
        nres: number of residues

    Returns:
        d: (nres+1, nres+1) int16 distance matrix, 1-based indexing
    """
    d = np.zeros((nres + 1, nres + 1), dtype=np.int16)
    for i in range(1, nres + 1):
        d[i, i] = 0
        for j in range(i + 1, nres + 1):
            dx = ca[0, i-1] - ca[0, j-1]
            dy = ca[1, i-1] - ca[1, j-1]
            dz = ca[2, i-1] - ca[2, j-1]
            dist = np.int16(_nint(10.0 * math.sqrt(dx*dx + dy*dy + dz*dz)))
            d[i, j] = dist
            d[j, i] = dist
    return d


def getp(x):
    """
    Metropolis acceptance probability.

    Args:
        x: alfa * dx (float32)

    Returns:
        probability in [0, 1] as float32 to match Fortran real precision
    """
    # Fortran getp: real function getp(x) — all float32
    x32 = np.float32(x)
    if x32 > np.float32(0.0):
        return np.float32(1.0)
    elif x32 < np.float32(-5.0):
        return np.float32(0.0)
    else:
        # Fortran exp() on real returns real (float32)
        return np.float32(np.exp(x32))


class RandomState:
    """
    Match gfortran intrinsic RAND() behavior.

    When compiled with gfortran, lean.f's `ran(seed)` resolves to the
    gfortran intrinsic RAND, NOT the user-defined ran.f LCG function.
    The intrinsic RAND(flag) with flag != 0 and flag != 1 reseeds and
    returns a fixed value. Since the common-block seed (124462287) never
    changes, every call returns the same constant: 0.0878167152.
    """

    def __init__(self, seed=124462287):
        self.seed = seed
        # gfortran RAND(124462287) always returns this constant
        self._constant = np.float32(8.78167152e-02)

    def ran(self):
        return self._constant
