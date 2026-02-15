use ndarray::{Array1, Array2, s};

/// Result of Kabsch superposition.
pub struct U3bResult {
    pub rms: f64,           // sum of weighted squared deviations (NOT root-mean-square)
    pub u: Array2<f64>,     // (3, 3) rotation matrix
    pub t: Array1<f64>,     // (3,) translation vector
    pub ier: i32,           // 0 = success, -1 = not unique, -2 = bad weights
}

/// Kabsch optimal rotation and translation.
///
/// Find rotation U and translation T such that U*X + T ~ Y,
/// minimizing sum(W * |U*X + T - Y|^2).
///
/// Uses SVD-based approach equivalent to the original Fortran eigenvalue method.
///
/// # Arguments
/// * `w` - (n,) weights
/// * `x` - (3, n) first coordinate set
/// * `y` - (3, n) second coordinate set
/// * `n` - number of atom pairs
/// * `mode` - 0 = RMS only, 1 = compute U, T
pub fn u3b(w: &[f64], x: &Array2<f64>, y: &Array2<f64>, n: usize, mode: i32) -> U3bResult {
    let mut u = Array2::eye(3);
    let mut t = Array1::zeros(3);
    let mut rms = 0.0;
    let mut ier: i32 = -1;

    if n < 1 {
        return U3bResult { rms, u, t, ier };
    }

    // Check weights
    let ww = &w[..n];
    if ww.iter().any(|&wi| wi < 0.0) {
        return U3bResult { rms, u, t, ier: -2 };
    }

    let wc: f64 = ww.iter().sum();
    if wc <= 0.0 {
        return U3bResult { rms, u, t, ier: -2 };
    }

    // Weighted centroids
    let mut xc = [0.0f64; 3];
    let mut yc = [0.0f64; 3];
    for i in 0..n {
        let wi = ww[i];
        for k in 0..3 {
            xc[k] += wi * x[[k, i]];
            yc[k] += wi * y[[k, i]];
        }
    }
    for k in 0..3 {
        xc[k] /= wc;
        yc[k] /= wc;
    }

    // E0 = sum of weighted squared distances from centroids
    let mut e0 = 0.0;
    for i in 0..n {
        let wi = ww[i];
        let mut sx2 = 0.0;
        let mut sy2 = 0.0;
        for k in 0..3 {
            let dx = x[[k, i]] - xc[k];
            let dy = y[[k, i]] - yc[k];
            sx2 += dx * dx;
            sy2 += dy * dy;
        }
        e0 += wi * (sx2 + sy2);
    }

    // Correlation matrix R = Y_centered * diag(W) * X_centered^T  (3x3)
    let mut r = Array2::zeros((3, 3));
    for i in 0..n {
        let wi = ww[i];
        for j in 0..3 {
            let dy = y[[j, i]] - yc[j];
            for k in 0..3 {
                let dx = x[[k, i]] - xc[k];
                r[[j, k]] += wi * dy * dx;
            }
        }
    }

    if mode == 0 {
        // RMS only via eigenvalues of R^T R
        let rtr = r.t().dot(&r);
        let eigenvalues = symmetric_eigenvalues_3x3(&rtr);
        let mut e: [f64; 3] = [0.0; 3];
        for i in 0..3 {
            e[i] = eigenvalues[i].max(0.0).sqrt();
        }
        // eigenvalues sorted descending
        e.sort_by(|a, b| b.partial_cmp(a).unwrap());

        let det_r = det3x3(&r);
        let d = if det_r >= 0.0 { e[2] } else { -e[2] };
        rms = e0 - 2.0 * (e[0] + e[1] + d);
        if rms < 0.0 {
            rms = 0.0;
        }
        ier = if e[1] <= e[0] * 1e-5 { -1 } else { 0 };

        // Translation from centroids (identity rotation)
        for k in 0..3 {
            t[k] = yc[k] - xc[k];
        }
        return U3bResult { rms, u, t, ier };
    }

    // SVD of R
    let (v_mat, s_vals, wt_mat) = svd_3x3(&r);

    // Handle reflection: ensure proper rotation
    let det_v = det3x3(&v_mat);
    let det_wt = det3x3(&wt_mat);
    let d = det_v * det_wt;

    let mut v = v_mat;
    let mut s = s_vals;
    if d < 0.0 {
        // Flip sign of last column of V and last singular value
        for k in 0..3 {
            v[[k, 2]] = -v[[k, 2]];
        }
        s[2] = -s[2];
    }

    // Rotation matrix U = V @ Wt
    u = v.dot(&wt_mat);

    // Translation
    let xc_arr = Array1::from_vec(xc.to_vec());
    let yc_arr = Array1::from_vec(yc.to_vec());
    t = &yc_arr - &u.dot(&xc_arr);

    // RMS = E0 - 2 * sum(singular values)
    rms = e0 - 2.0 * s.iter().sum::<f64>();
    if rms < 0.0 {
        rms = 0.0;
    }

    ier = if s[1] <= s[0].abs() * 1e-5 { -1 } else { 0 };

    U3bResult { rms, u, t, ier }
}

/// Apply rotation and translation in-place: x = U @ x + T
pub fn transrotate(x: &mut Array2<f64>, u: &Array2<f64>, t: &Array1<f64>) {
    let n = x.ncols();
    for i in 0..n {
        let col = x.slice(s![.., i]).to_owned();
        let rotated = u.dot(&col) + t;
        x[[0, i]] = rotated[0];
        x[[1, i]] = rotated[1];
        x[[2, i]] = rotated[2];
    }
}

/// Compute rotation matrix and translation vector from alignment blocks.
///
/// Extracts aligned CA pairs from the alignment blocks, computes Kabsch
/// superposition (U, T) such that U*ca1 + T ~ ca2 for the aligned positions.
///
/// Returns None if no aligned pairs exist.
pub fn compute_transform(
    ca1: &Array2<f64>,  // (3, nres1)
    ca2: &Array2<f64>,  // (3, nres2)
    blocks: &[crate::AlignmentBlock],
) -> Option<(Array2<f64>, Array1<f64>)> {
    // Count total aligned pairs
    let n: usize = blocks.iter().map(|b| (b.r1 - b.l1 + 1) as usize).sum();
    if n == 0 {
        return None;
    }

    // Extract aligned coordinates
    let mut x = Array2::zeros((3, n));
    let mut y = Array2::zeros((3, n));
    let mut idx = 0;
    for b in blocks {
        let len = (b.r1 - b.l1 + 1) as usize;
        for k in 0..len {
            let i1 = (b.l1 as usize - 1) + k; // 0-based
            let i2 = (b.l2 as usize - 1) + k;
            for d in 0..3 {
                x[[d, idx]] = ca1[[d, i1]];
                y[[d, idx]] = ca2[[d, i2]];
            }
            idx += 1;
        }
    }

    let w: Vec<f64> = vec![1.0; n];
    let result = u3b(&w, &x, &y, n, 1);
    Some((result.u, result.t))
}

// --- Internal linear algebra for 3x3 matrices ---

fn det3x3(m: &Array2<f64>) -> f64 {
    m[[0, 0]] * (m[[1, 1]] * m[[2, 2]] - m[[1, 2]] * m[[2, 1]])
        - m[[0, 1]] * (m[[1, 0]] * m[[2, 2]] - m[[1, 2]] * m[[2, 0]])
        + m[[0, 2]] * (m[[1, 0]] * m[[2, 1]] - m[[1, 1]] * m[[2, 0]])
}

/// Compute eigenvalues of a 3x3 symmetric matrix using Cardano's formula.
fn symmetric_eigenvalues_3x3(m: &Array2<f64>) -> [f64; 3] {
    let a = m[[0, 0]];
    let b = m[[1, 1]];
    let c = m[[2, 2]];
    let d = m[[0, 1]];
    let e = m[[0, 2]];
    let f = m[[1, 2]];

    // Characteristic polynomial: -lambda^3 + p*lambda^2 + q*lambda + r = 0
    // trace
    let p = a + b + c;
    // sum of 2x2 minors
    let q = d * d + e * e + f * f - a * b - a * c - b * c;
    // determinant
    let r = a * b * c + 2.0 * d * e * f - a * f * f - b * e * e - c * d * d;

    // Solve depressed cubic: t^3 + pt + q = 0 where we substitute lambda = t + p/3
    let p3 = p / 3.0;
    let qq = q + p * p / 3.0;  // coefficient of depressed cubic (negated from standard)
    let rr = r + p * qq / 3.0 - 2.0 * p3 * p3 * p3 + p3 * q;

    // Use trigonometric solution for 3 real roots
    let qq3 = qq / 3.0;
    let disc = qq3 * qq3 * qq3 + (rr / 2.0) * (rr / 2.0);

    if disc <= 0.0 {
        // Three real roots
        let mag = (-qq3).max(0.0).sqrt();
        let theta = if mag > 0.0 {
            ((-rr / 2.0) / (mag * mag * mag)).clamp(-1.0, 1.0).acos() / 3.0
        } else {
            0.0
        };

        let two_pi_3 = 2.0 * std::f64::consts::PI / 3.0;
        let mut eigs = [
            2.0 * mag * theta.cos() + p3,
            2.0 * mag * (theta - two_pi_3).cos() + p3,
            2.0 * mag * (theta + two_pi_3).cos() + p3,
        ];
        eigs.sort_by(|a, b| b.partial_cmp(a).unwrap());
        eigs
    } else {
        // One real root, two complex (shouldn't happen for symmetric matrix)
        let sqrt_disc = disc.sqrt();
        let s = cbrt(-rr / 2.0 + sqrt_disc);
        let t = cbrt(-rr / 2.0 - sqrt_disc);
        let root = s + t + p3;
        [root, root, root]
    }
}

fn cbrt(x: f64) -> f64 {
    if x >= 0.0 {
        x.powf(1.0 / 3.0)
    } else {
        -(-x).powf(1.0 / 3.0)
    }
}

/// SVD of a 3x3 matrix using iterative Jacobi method.
///
/// Returns (U, S, Vt) where M = U * diag(S) * Vt
/// S is sorted descending by absolute value.
fn svd_3x3(m: &Array2<f64>) -> (Array2<f64>, [f64; 3], Array2<f64>) {
    // Compute M^T M
    let mtm = m.t().dot(m);

    // Eigendecomposition of M^T M using Jacobi iteration
    let (eigenvalues, v) = jacobi_eigen_3x3(&mtm);

    // Singular values = sqrt of eigenvalues (sorted descending)
    let mut indices = [0usize, 1, 2];
    indices.sort_by(|&a, &b| eigenvalues[b].partial_cmp(&eigenvalues[a]).unwrap());

    let mut s = [0.0f64; 3];
    let mut v_sorted = Array2::zeros((3, 3));
    for (new_i, &old_i) in indices.iter().enumerate() {
        s[new_i] = eigenvalues[old_i].max(0.0).sqrt();
        for k in 0..3 {
            v_sorted[[k, new_i]] = v[[k, old_i]];
        }
    }

    // U = M * V * S^{-1}
    let mut u_mat = Array2::zeros((3, 3));
    for i in 0..3 {
        if s[i] > 1e-15 {
            let v_col = v_sorted.column(i).to_owned();
            let mv = m.dot(&v_col);
            for k in 0..3 {
                u_mat[[k, i]] = mv[k] / s[i];
            }
        }
    }

    // Handle degenerate case: if s[2] ~ 0, compute u[:,2] as cross product
    if s[2] < 1e-15 {
        let u0 = [u_mat[[0, 0]], u_mat[[1, 0]], u_mat[[2, 0]]];
        let u1 = [u_mat[[0, 1]], u_mat[[1, 1]], u_mat[[2, 1]]];
        u_mat[[0, 2]] = u0[1] * u1[2] - u0[2] * u1[1];
        u_mat[[1, 2]] = u0[2] * u1[0] - u0[0] * u1[2];
        u_mat[[2, 2]] = u0[0] * u1[1] - u0[1] * u1[0];
    }

    let vt = v_sorted.t().to_owned();
    (u_mat, s, vt)
}

/// Jacobi eigendecomposition of a 3x3 symmetric matrix.
/// Returns (eigenvalues, eigenvector_matrix).
fn jacobi_eigen_3x3(m: &Array2<f64>) -> ([f64; 3], Array2<f64>) {
    let mut a = m.clone();
    let mut v = Array2::eye(3);
    let max_iter = 100;

    for _ in 0..max_iter {
        // Find largest off-diagonal element
        let mut max_val = 0.0f64;
        let mut p = 0usize;
        let mut q = 1usize;
        for i in 0..3 {
            for j in (i + 1)..3 {
                if a[[i, j]].abs() > max_val {
                    max_val = a[[i, j]].abs();
                    p = i;
                    q = j;
                }
            }
        }

        if max_val < 1e-15 {
            break;
        }

        // Compute rotation
        let theta = if (a[[q, q]] - a[[p, p]]).abs() < 1e-30 {
            std::f64::consts::PI / 4.0
        } else {
            0.5 * (2.0 * a[[p, q]] / (a[[p, p]] - a[[q, q]])).atan()
        };

        let c = theta.cos();
        let s = theta.sin();

        // Apply Jacobi rotation: A' = G^T A G
        let mut a_new = a.clone();

        // Update diagonal elements
        a_new[[p, p]] = c * c * a[[p, p]] + 2.0 * s * c * a[[p, q]] + s * s * a[[q, q]];
        a_new[[q, q]] = s * s * a[[p, p]] - 2.0 * s * c * a[[p, q]] + c * c * a[[q, q]];
        a_new[[p, q]] = 0.0;
        a_new[[q, p]] = 0.0;

        // Update off-diagonal elements
        for r in 0..3 {
            if r != p && r != q {
                let arp = c * a[[r, p]] + s * a[[r, q]];
                let arq = -s * a[[r, p]] + c * a[[r, q]];
                a_new[[r, p]] = arp;
                a_new[[p, r]] = arp;
                a_new[[r, q]] = arq;
                a_new[[q, r]] = arq;
            }
        }

        a = a_new;

        // Update eigenvector matrix
        let mut v_new = v.clone();
        for r in 0..3 {
            v_new[[r, p]] = c * v[[r, p]] + s * v[[r, q]];
            v_new[[r, q]] = -s * v[[r, p]] + c * v[[r, q]];
        }
        v = v_new;
    }

    ([a[[0, 0]], a[[1, 1]], a[[2, 2]]], v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn test_u3b_identity() {
        // Non-collinear points so SVD is well-determined
        let w = vec![1.0, 1.0, 1.0, 1.0];
        let coords = array![
            [0.0, 1.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 1.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let result = u3b(&w, &coords, &coords, 4, 1);
        assert!(result.rms < 1e-10, "rms should be ~0, got {}", result.rms);
        assert_eq!(result.ier, 0);
    }

    #[test]
    fn test_u3b_pure_translation() {
        let w = vec![1.0, 1.0, 1.0, 1.0];
        let x = array![
            [0.0, 1.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 1.0],
            [0.0, 0.0, 0.0, 0.0],
        ];
        let y = &x + 5.0;
        let result = u3b(&w, &x, &y, 4, 1);
        assert!(result.rms < 1e-10);
        // Translation should be ~[5,5,5]
        for k in 0..3 {
            assert!((result.t[k] - 5.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_transrotate() {
        let u = Array2::eye(3);
        let t = array![1.0, 2.0, 3.0];
        let mut x = array![
            [0.0, 1.0],
            [0.0, 1.0],
            [0.0, 1.0],
        ];
        transrotate(&mut x, &u, &t);
        assert!((x[[0, 0]] - 1.0).abs() < 1e-10);
        assert!((x[[1, 0]] - 2.0).abs() < 1e-10);
        assert!((x[[2, 0]] - 3.0).abs() < 1e-10);
    }
}
