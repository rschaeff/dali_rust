use ndarray::Array2;

/// Build similarity score matrix for NW alignment.
///
/// score[i,j] = max(0, nint((rcut - dist(x[:,i], y[:,j])) * 10))
pub fn filltable_maxsim(x: &Array2<f64>, y: &Array2<f64>, rcut: f64) -> Array2<i16> {
    let nx = x.ncols();
    let ny = y.ncols();
    let mut table = Array2::zeros((nx, ny));

    for i in 0..nx {
        for j in 0..ny {
            let dx = x[[0, i]] - y[[0, j]];
            let dy = x[[1, i]] - y[[1, j]];
            let dz = x[[2, i]] - y[[2, j]];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < rcut {
                // Use Rust's f64::round() which matches rint() for this usage
                let s = ((rcut - r) * 10.0).round() as i16;
                if s > 0 {
                    table[[i, j]] = s;
                }
            }
        }
    }

    table
}

/// Anti-diagonal Needleman-Wunsch alignment maximizing similarity.
/// Zero gap penalty.
///
/// Translated from subfitzfast.f nw_maxsim().
///
/// Returns (score, ali) where ali[i] = j (0-based) or -1 for unaligned.
pub fn nw_maxsim(m: usize, n: usize, table0: &Array2<i16>) -> (i32, Vec<i32>) {
    if m < 1 || n < 1 {
        return (0, vec![-1; m]);
    }

    // Create 1-based padded table
    let mut t = Array2::<i16>::zeros((m + 1, n + 1));
    for i in 0..m {
        for j in 0..n {
            t[[i + 1, j + 1]] = table0[[i, j]];
        }
    }

    let mut trace = Array2::<i16>::zeros((m + 1, n + 1));

    // Initialize boundary entries
    trace[[1, 1]] = 0;
    if m >= 2 {
        trace[[2, 1]] = 1;
    }
    if n >= 2 {
        trace[[1, 2]] = 2;
    }

    let max_k = m.max(n) + 1;
    let mut shift = vec![vec![0i64; max_k]; 3];
    let mut diag0: usize = 2;
    let mut diag1: usize = 1;
    let mut diag2: usize = 0;

    shift[diag2][0] = t[[1, 1]] as i64;
    if m >= 2 {
        shift[diag1][0] = (t[[2, 1]] as i64).max(t[[1, 1]] as i64);
    }
    if n >= 2 {
        shift[diag1][1] = (t[[1, 2]] as i64).max(t[[1, 1]] as i64);
    }

    // Main anti-diagonal fill
    for y0 in 3..=(m + n - 1) {
        let k0 = if y0 > n { y0 - n } else { 0 };
        let k1 = (m - 1).min(y0 - 1);

        if k0 == k1 {
            // Single cell on this anti-diagonal — check all available options
            let x = 1 + k0;
            let y = y0 - k0;
            if x <= m && y >= 1 && y <= n {
                let s = t[[x, y]] as i64;
                let b = shift[diag1][k0] + s;              // gap in x (from left)
                let a = if k0 > 0 { shift[diag1][k0 - 1] + s } else { i64::MIN }; // gap in y
                let c = if k0 > 0 { shift[diag2][k0 - 1] + s } else { i64::MIN }; // diagonal

                if c >= a && c >= b {
                    trace[[x, y]] = 3;
                    shift[diag0][k0] = c;
                } else if b >= a {
                    trace[[x, y]] = 2;
                    shift[diag0][k0] = b;
                } else {
                    trace[[x, y]] = 1;
                    shift[diag0][k0] = a;
                }
            }
        } else {
            // Left boundary
            let x = 1 + k0;
            let y = y0 - k0;
            if x <= m && y >= 1 && y <= n {
                trace[[x, y]] = 2;
                shift[diag0][k0] = shift[diag1][k0] + t[[x, y]] as i64;
            }

            // Right boundary
            let x = 1 + k1;
            let y = y0 - k1;
            if x <= m && y >= 1 && y <= n && k1 > 0 {
                trace[[x, y]] = 1;
                shift[diag0][k1] = shift[diag1][k1 - 1] + t[[x, y]] as i64;
            }

            // Interior cells
            for k in (k0 + 1)..k1 {
                let x = 1 + k;
                let y = y0 - k;
                if x > m || y < 1 || y > n {
                    continue;
                }
                let s = t[[x, y]] as i64;
                let a = shift[diag1][k - 1]; // gap in y
                let b = shift[diag1][k];     // gap in x
                let c = shift[diag2][k - 1]; // diagonal

                if c >= a && c >= b {
                    trace[[x, y]] = 3;
                    shift[diag0][k] = c + s;
                } else if b >= a && b >= c {
                    trace[[x, y]] = 2;
                    shift[diag0][k] = b;
                } else {
                    trace[[x, y]] = 1;
                    shift[diag0][k] = a;
                }
            }
        }

        // Roll diagonal indices
        diag0 = (diag0 + 1) % 3;
        diag1 = (diag1 + 1) % 3;
        diag2 = (diag2 + 1) % 3;
    }

    // Backtrack alignment (1-based)
    let mut ali_1 = vec![0i32; m + 1];
    let mut x = m;
    let mut y = n;
    while x > 1 || y > 1 {
        let tr = trace[[x, y]];
        match tr {
            1 => x -= 1,
            2 => y -= 1,
            3 => {
                if t[[x, y]] > 0 {
                    ali_1[x] = y as i32;
                }
                x -= 1;
                y -= 1;
            }
            _ => {
                x = 0;
                y = 0;
            }
        }
    }

    // Compute score
    let mut score: i32 = 0;
    for i in 1..=m {
        if ali_1[i] > 0 {
            score += t[[i, ali_1[i] as usize]] as i32;
        }
    }

    // Convert to 0-based output
    let mut ali = vec![-1i32; m];
    for i in 1..=m {
        if ali_1[i] > 0 {
            ali[i - 1] = ali_1[i] - 1;
        }
    }

    (score, ali)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn test_filltable_maxsim_identical() {
        let x = array![[0.0, 1.0], [0.0, 0.0], [0.0, 0.0]];
        let table = filltable_maxsim(&x, &x, 4.0);
        // Self-distance is 0, score = nint((4.0 - 0.0) * 10) = 40
        assert_eq!(table[[0, 0]], 40);
        assert_eq!(table[[1, 1]], 40);
        // Cross: dist = 1.0, score = nint((4.0 - 1.0) * 10) = 30
        assert_eq!(table[[0, 1]], 30);
        assert_eq!(table[[1, 0]], 30);
    }

    #[test]
    fn test_nw_maxsim_identity() {
        // 3x3 identity-like table: diagonal = 40, off-diagonal = 0.
        // Fortran/Python behavior: trace[1,1]=0, so cell (1,1) is never aligned
        // during backtracking. Only cells (2,2) and (3,3) get aligned.
        let table = array![
            [40i16, 0,  0],
            [0,    40,  0],
            [0,     0, 40],
        ];
        let (score, ali) = nw_maxsim(3, 3, &table);
        assert_eq!(ali[0], -1); // cell (1,1) not aligned (Fortran behavior)
        assert_eq!(ali[1], 1);
        assert_eq!(ali[2], 2);
        assert_eq!(score, 80);
    }
}
