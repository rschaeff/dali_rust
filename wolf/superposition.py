"""
Kabsch superposition (u3b) and coordinate transformation.

u3b translated from u3b-8.f — Kabsch (1976, 1978).
"""

import numpy as np


def u3b(w, x, y, n, mode=1):
    """
    Kabsch optimal rotation and translation.

    Find rotation U and translation T such that U*X + T ≈ Y,
    minimizing sum(W * |U*X + T - Y|^2).

    Uses SVD-based approach equivalent to the original Fortran eigenvalue method.

    Args:
        w: (n,) weights
        x: (3, n) first coordinate set
        y: (3, n) second coordinate set
        n: number of atom pairs
        mode: 0 = RMS only, 1 = compute U, T

    Returns:
        rms: sum of weighted squared deviations (NOT root-mean-square)
        u: (3, 3) rotation matrix
        t: (3,) translation vector
        ier: error code (0 = success, -1 = not unique, -2 = bad weights)
    """
    u = np.eye(3, dtype=np.float64)
    t = np.zeros(3, dtype=np.float64)
    rms = 0.0
    ier = -1

    if n < 1:
        return rms, u, t, ier

    # Check weights
    ww = w[:n].astype(np.float64)
    if np.any(ww < 0):
        ier = -2
        return rms, u, t, ier

    wc = np.sum(ww)
    if wc <= 0:
        ier = -2
        return rms, u, t, ier

    # Weighted centroids
    xc = np.sum(ww[None, :] * x[:, :n], axis=1) / wc
    yc = np.sum(ww[None, :] * y[:, :n], axis=1) / wc

    # Centered coordinates
    xm = x[:, :n] - xc[:, None]  # (3, n)
    ym = y[:, :n] - yc[:, None]  # (3, n)

    # E0 = sum of weighted squared distances from centroids
    e0 = np.sum(ww * np.sum(xm ** 2, axis=0)) + np.sum(ww * np.sum(ym ** 2, axis=0))

    # Correlation matrix R = Y * diag(W) * X^T
    R = (ym * ww[None, :]) @ xm.T  # (3, 3)

    # Determinant of R
    det_R = np.linalg.det(R)
    sigma = det_R

    if mode == 0:
        # RMS only via eigenvalues of R^T R
        RtR = R.T @ R
        eigenvalues = np.linalg.eigvalsh(RtR)
        eigenvalues = np.sort(eigenvalues)[::-1]  # descending
        e = np.sqrt(np.maximum(eigenvalues, 0))
        d = e[2] if sigma >= 0 else -e[2]
        rms = e0 - 2.0 * (e[0] + e[1] + d)
        if rms < 0:
            rms = 0.0
        ier = 0
        if e[1] <= e[0] * 1e-5:
            ier = -1
        # Translation from centroids (with identity rotation)
        t = yc - xc
        return rms, u, t, ier

    # SVD of R
    V, S, Wt = np.linalg.svd(R)

    # Handle reflection: ensure proper rotation
    d = np.linalg.det(V) * np.linalg.det(Wt)
    if d < 0:
        # Flip the sign of the last singular vector
        V[:, -1] *= -1
        S[-1] *= -1

    # Rotation matrix
    u = V @ Wt

    # Translation
    t = yc - u @ xc

    # RMS = E0 - 2 * sum(singular values, with sign correction)
    rms = e0 - 2.0 * np.sum(S)
    if rms < 0:
        rms = 0.0

    ier = 0
    if len(S) >= 2 and S[1] <= S[0] * 1e-5:
        ier = -1

    return rms, u, t, ier


def transrotate(x, u, t):
    """
    Apply rotation and translation in-place: x = U @ x + T.

    Args:
        x: (3, n) coordinates, modified in-place
        u: (3, 3) rotation matrix
        t: (3,) translation vector
    """
    n = x.shape[1]
    for i in range(n):
        tmp = t + u @ x[:, i]
        x[:, i] = tmp
