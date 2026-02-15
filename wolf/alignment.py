"""
Needleman-Wunsch dynamic programming and fitz iterative refinement.

Translated from subfitzfast.f.
"""

import numpy as np
from .superposition import u3b, transrotate


def filltable_maxsim(x, y, rcut):
    """
    Build similarity score matrix for NW alignment.

    score[i,j] = max(0, nint((rcut - dist(x[:,i], y[:,j])) * 10))

    Args:
        x: (3, nx) first coordinate set
        y: (3, ny) second coordinate set
        rcut: distance cutoff

    Returns:
        table0: (nx, ny) integer score matrix (1-based internally, returned 0-based)
    """
    nx = x.shape[1]
    ny = y.shape[1]
    table0 = np.zeros((nx, ny), dtype=np.int16)

    for i in range(nx):
        for j in range(ny):
            dx = x[0, i] - y[0, j]
            dy = x[1, i] - y[1, j]
            dz = x[2, i] - y[2, j]
            r = np.sqrt(dx * dx + dy * dy + dz * dz)
            if r < rcut:
                table0[i, j] = int(np.rint((rcut - r) * 10.0))
            # else: stays 0

    return table0


def nw_maxsim(m, n, table0):
    """
    Anti-diagonal Needleman-Wunsch alignment maximizing similarity.
    Zero gap penalty.

    Translated from subfitzfast.f nw_maxsim().
    Uses 1-based indexing internally for exact Fortran correspondence.

    Args:
        m: number of residues in sequence 1
        n: number of residues in sequence 2
        table0: (m, n) score matrix (0-based)

    Returns:
        score: total alignment score
        ali: (m,) alignment array. ali[i] = j (0-based) or -1 for unaligned.
    """
    if m < 1 or n < 1:
        return 0, np.full(m, -1, dtype=np.int32)

    # Create 1-based padded table
    t = np.zeros((m + 1, n + 1), dtype=np.int16)
    t[1:m + 1, 1:n + 1] = table0

    trace = np.zeros((m + 1, n + 1), dtype=np.int16)

    # Initialize boundary entries
    trace[1, 1] = 0
    if m >= 2:
        trace[2, 1] = 1
    if n >= 2:
        trace[1, 2] = 2

    max_k = max(m, n) + 1
    shift = np.zeros((3, max_k), dtype=np.int64)
    diag0, diag1, diag2 = 2, 1, 0

    shift[diag2, 0] = int(t[1, 1])
    if m >= 2:
        shift[diag1, 0] = max(int(t[2, 1]), int(t[1, 1]))
    if n >= 2:
        shift[diag1, 1] = max(int(t[1, 2]), int(t[1, 1]))

    # Main anti-diagonal fill
    for y0 in range(3, m + n):
        k0 = max(0, y0 - n)     # boundary: x = 1
        k1 = min(m - 1, y0 - 1)  # boundary: y = 1

        # Left boundary (x minimal, y maximal)
        x = 1 + k0
        y = y0 - k0
        trace[x, y] = 2
        shift[diag0, k0] = int(shift[diag1, k0]) + int(t[x, y])

        # Right boundary (x maximal, y minimal)
        x = 1 + k1
        y = y0 - k1
        trace[x, y] = 1
        shift[diag0, k1] = int(shift[diag1, k1 - 1]) + int(t[x, y])

        # Interior cells
        for k in range(k0 + 1, k1):
            x = 1 + k
            y = y0 - k
            s = int(t[x, y])
            a = int(shift[diag1, k - 1])   # gap in y (x advances)
            b = int(shift[diag1, k])         # gap in x (y advances)
            c = int(shift[diag2, k - 1])     # diagonal (match)

            if c >= a and c >= b:
                trace[x, y] = 3
                shift[diag0, k] = c + s
            elif b >= a and b >= c:
                trace[x, y] = 2
                shift[diag0, k] = b
            else:
                trace[x, y] = 1
                shift[diag0, k] = a

        # Roll diagonal indices
        diag0 = (diag0 + 1) % 3
        diag1 = (diag1 + 1) % 3
        diag2 = (diag2 + 1) % 3

    # Backtrack alignment (1-based)
    ali_1 = np.zeros(m + 1, dtype=np.int32)
    x, y = m, n
    while x > 1 or y > 1:
        tr = trace[x, y]
        if tr == 1:
            x -= 1
        elif tr == 2:
            y -= 1
        elif tr == 3:
            if t[x, y] > 0:
                ali_1[x] = y
            x -= 1
            y -= 1
        else:
            x = 0
            y = 0

    # Compute score
    score = 0
    for i in range(1, m + 1):
        if ali_1[i] > 0:
            score += int(t[i, ali_1[i]])

    # Convert to 0-based output
    ali = np.full(m, -1, dtype=np.int32)
    for i in range(1, m + 1):
        if ali_1[i] > 0:
            ali[i - 1] = ali_1[i] - 1

    return score, ali


def getut(nx, ali, x, y):
    """
    Compute optimal rotation/translation from alignment.

    Args:
        nx: length of x
        ali: (nx,) alignment, ali[i] = j (0-based) or -1
        x: (3, nx)
        y: (3, ny)

    Returns:
        u: (3, 3) rotation matrix
        t: (3,) translation vector
        lali: alignment length
        rms: RMSD
    """
    # Extract aligned pairs
    pairs = [(i, ali[i]) for i in range(nx) if ali[i] >= 0]
    lali = len(pairs)

    if lali == 0:
        return np.eye(3), np.zeros(3), 0, 0.0

    ux = np.zeros((3, lali))
    uy = np.zeros((3, lali))
    w = np.ones(lali)

    for k, (i, j) in enumerate(pairs):
        ux[:, k] = x[:, i]
        uy[:, k] = y[:, j]

    ssq, u, t, ier = u3b(w, ux, uy, lali, mode=1)
    rms_val = np.sqrt(ssq / max(1, lali))

    return u, t, lali, rms_val


def fitz(x, y, rcut, maxiter):
    """
    Iterative superposition refinement.

    Repeatedly: score matrix → NW alignment → Kabsch superposition,
    until alignment freezes or maxiter reached.

    Modifies x in-place (applies rotation/translation each iteration).

    Args:
        x: (3, nx) first coordinate set (modified in-place)
        y: (3, ny) second coordinate set
        rcut: distance cutoff for scoring
        maxiter: maximum iterations

    Returns:
        ali: (nx,) final alignment, ali[i] = j (0-based) or -1
        rms: final RMSD
        lali: alignment length
        niter: number of iterations performed
    """
    nx = x.shape[1]
    ny = y.shape[1]

    if nx < 1 or ny < 1:
        return np.full(max(nx, 0), -1, dtype=np.int32), 0.0, 0, 0

    ali = np.full(nx, -1, dtype=np.int32)
    rms = 0.0
    lali = 0
    niter = 0
    frozen = False

    while niter < maxiter and not frozen:
        niter += 1

        # Build score table
        table0 = filltable_maxsim(x, y, rcut)

        # Save old alignment
        ali_old = ali.copy()

        # NW alignment
        score, ali = nw_maxsim(nx, ny, table0)

        # Check frozen
        frozen = np.array_equal(ali, ali_old)

        # Superpose
        u, t, lali, rms = getut(nx, ali, x, y)
        transrotate(x, u, t)

    return ali, rms, lali, niter
