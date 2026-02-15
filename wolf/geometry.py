"""
SSE geometry: vec(), twist(), preparex().

Translated from wolf_original.f.
"""

import numpy as np
import math


def vec(ca, nres, nseg, segmentrange):
    """
    Compute SSE midpoints and direction vectors from CA coordinates.

    For each SSE, split residues into N-half and C-half:
      mid = average of (N-half center + C-half center) / 2
      dir = C-half center - N-half center

    Args:
        ca: (3, nres) CA coordinates (0-indexed)
        nres: number of residues
        nseg: number of SSEs
        segmentrange: (2, nseg) start/end residue indices (1-based)

    Returns:
        midpoint: (3, nseg) SSE midpoints
        direction: (3, nseg) SSE direction vectors
    """
    midpoint = np.zeros((3, nseg), dtype=np.float64)
    direction = np.zeros((3, nseg), dtype=np.float64)

    for iseg in range(nseg):
        left = segmentrange[0, iseg]   # 1-based
        rite = segmentrange[1, iseg]   # 1-based
        mid_res = (left + rite) // 2   # Fortran integer division
        l = mid_res - left + 1         # N-half count
        r = rite - mid_res + 1         # C-half count (note: mid_res is shared)

        # N-half: residues left..mid_res (inclusive, 1-based)
        nmid = np.zeros(3)
        for ires in range(left, mid_res + 1):  # 1-based inclusive
            nmid += ca[:, ires - 1] / l

        # C-half: residues mid_res..rite (inclusive, 1-based)
        cmid = np.zeros(3)
        for ires in range(mid_res, rite + 1):  # 1-based inclusive
            cmid += ca[:, ires - 1] / r

        midpoint[:, iseg] = (nmid + cmid) / 2.0
        direction[:, iseg] = cmid - nmid

    return midpoint, direction


def twist(x, looplen):
    """
    Transform coordinates to canonical frame.

    Puts x[:,0] at origin, x[:,1] along +y axis, x[:,2] in positive yz-plane.
    Modifies x in-place.

    Args:
        x: (3, looplen) coordinate array
        looplen: number of points
    """
    PI = 3.141592653589793

    # Set origin: subtract x[:,0] from all points
    origin = x[:, 0].copy()
    for i in range(looplen - 1, -1, -1):
        x[:, i] -= origin

    # Rotate around x-axis to put x[:,1] in xy-plane
    # If y ~= 0, first rotate 90° around z
    if abs(x[1, 1]) < 1e-6:
        for i in range(1, looplen):
            y0 = -x[1, i]
            y1 = x[0, i]
            y2 = x[2, i]
            x[0, i] = y0
            x[1, i] = y1
            x[2, i] = y2

    # If y ~= 0 and x ~= 0, rotate 90° around x
    if abs(x[1, 1]) < 1e-6 and abs(x[0, 1]) < 1e-6:
        for i in range(1, looplen):
            y1 = -x[2, i]
            y2 = x[1, i]
            y0 = x[0, i]
            x[0, i] = y0
            x[1, i] = y1
            x[2, i] = y2

    # Rotate around x-axis
    if abs(x[1, 1]) > 1e-6:
        u = math.atan(x[2, 1] / x[1, 1])
    else:
        u = 0.0
    sinu = math.sin(u)
    cosu = math.cos(u)
    for i in range(1, looplen):
        y1 = cosu * x[1, i] + sinu * x[2, i]
        y2 = -sinu * x[1, i] + cosu * x[2, i]
        y0 = x[0, i]
        x[0, i] = y0
        x[1, i] = y1
        x[2, i] = y2

    # Rotate around z-axis
    if abs(x[1, 1]) > 1e-6:
        u = -math.atan(x[0, 1] / x[1, 1])
    else:
        u = 0.0
    if x[1, 1] < 0.0:
        u = PI + u
    sinu = math.sin(u)
    cosu = math.cos(u)
    for i in range(1, looplen):
        y0 = cosu * x[0, i] + sinu * x[1, i]
        y1 = -sinu * x[0, i] + cosu * x[1, i]
        y2 = x[2, i]
        x[0, i] = y0
        x[1, i] = y1
        x[2, i] = y2

    # Rotate around y-axis to put x[:,2] in yz-plane
    # If z ~= 0, first rotate 90° around y
    if abs(x[2, 2]) < 1e-6:
        for i in range(1, looplen):
            y0 = x[2, i]
            y1 = x[1, i]
            y2 = -x[0, i]
            x[0, i] = y0
            x[1, i] = y1
            x[2, i] = y2

    if x[2, 2] != 0.0:
        u = math.atan(x[0, 2] / x[2, 2])
    else:
        u = 0.0
    sinu = math.sin(u)
    cosu = math.cos(u)
    for i in range(2, looplen):
        y2 = cosu * x[2, i] + sinu * x[0, i]
        y0 = -sinu * x[2, i] + cosu * x[0, i]
        y1 = x[1, i]
        x[0, i] = y0
        x[1, i] = y1
        x[2, i] = y2

    # If z < 0, rotate 180° around y
    if x[2, 2] < 0.0:
        for i in range(1, looplen):
            x[0, i] = -x[0, i]
            x[2, i] = -x[2, i]


def preparex(iseg, jseg, nseg, midpoint, direction):
    """
    Build reference frame array for twist().

    Layout:
      x[:,0] = midpoint[:,iseg]
      x[:,1] = midpoint[:,iseg] + direction[:,iseg]
      x[:,2] = midpoint[:,jseg]
      x[:,3:3+nseg] = midpoint[:,k] - direction[:,k]  for k=0..nseg-1
      x[:,3+nseg:3+2*nseg] = midpoint[:,k] + direction[:,k]  for k=0..nseg-1

    Args:
        iseg, jseg: SSE indices (0-based)
        nseg: number of SSEs
        midpoint: (3, nseg)
        direction: (3, nseg)

    Returns:
        x: (3, 3 + 2*nseg) coordinate array
    """
    n = 3 + 2 * nseg
    x = np.zeros((3, n), dtype=np.float64)

    x[:, 0] = midpoint[:, iseg]
    x[:, 1] = midpoint[:, iseg] + direction[:, iseg]
    x[:, 2] = midpoint[:, jseg]

    for kseg in range(nseg):
        x[:, 3 + kseg] = midpoint[:, kseg] - direction[:, kseg]
        x[:, 3 + nseg + kseg] = midpoint[:, kseg] + direction[:, kseg]

    return x


def compute_neidist(nseg, midpoint):
    """
    Compute pairwise Euclidean distances between SSE midpoints.

    Args:
        nseg: number of SSEs
        midpoint: (3, nseg) SSE midpoints

    Returns:
        neidist: (nseg, nseg) distance matrix
    """
    neidist = np.zeros((nseg, nseg), dtype=np.float64)
    for i in range(nseg):
        for j in range(nseg):
            if i != j:
                d = midpoint[:, i] - midpoint[:, j]
                neidist[i, j] = np.sqrt(np.sum(d * d))
    return neidist
