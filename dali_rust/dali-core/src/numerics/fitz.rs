use ndarray::Array2;

use super::kabsch::{u3b, transrotate};
use super::nw::{filltable_maxsim, nw_maxsim};

/// Result of iterative superposition refinement.
pub struct FitzResult {
    pub ali: Vec<i32>,     // alignment, ali[i] = j (0-based) or -1
    pub rms: f64,          // final RMSD
    pub lali: usize,       // alignment length
    pub niter: usize,      // iterations performed
}

/// Compute optimal rotation/translation from alignment.
///
/// Returns (u, t, lali, rms)
fn getut(
    nx: usize,
    ali: &[i32],
    x: &Array2<f64>,
    y: &Array2<f64>,
) -> (ndarray::Array2<f64>, ndarray::Array1<f64>, usize, f64) {
    // Extract aligned pairs
    let pairs: Vec<(usize, usize)> = (0..nx)
        .filter(|&i| ali[i] >= 0)
        .map(|i| (i, ali[i] as usize))
        .collect();
    let lali = pairs.len();

    if lali == 0 {
        return (
            Array2::eye(3),
            ndarray::Array1::zeros(3),
            0,
            0.0,
        );
    }

    let mut ux = Array2::zeros((3, lali));
    let mut uy = Array2::zeros((3, lali));
    let w = vec![1.0; lali];

    for (k, &(i, j)) in pairs.iter().enumerate() {
        for dim in 0..3 {
            ux[[dim, k]] = x[[dim, i]];
            uy[[dim, k]] = y[[dim, j]];
        }
    }

    let result = u3b(&w, &ux, &uy, lali, 1);
    let rms_val = (result.rms / lali.max(1) as f64).sqrt();

    (result.u, result.t, lali, rms_val)
}

/// Iterative superposition refinement.
///
/// Repeatedly: score matrix -> NW alignment -> Kabsch superposition,
/// until alignment freezes or maxiter reached.
///
/// Modifies x in-place (applies rotation/translation each iteration).
pub fn fitz(x: &mut Array2<f64>, y: &Array2<f64>, rcut: f64, maxiter: usize) -> FitzResult {
    let nx = x.ncols();
    let ny = y.ncols();

    if nx < 1 || ny < 1 {
        return FitzResult {
            ali: vec![-1; nx],
            rms: 0.0,
            lali: 0,
            niter: 0,
        };
    }

    let mut ali = vec![-1i32; nx];
    let mut rms = 0.0;
    let mut lali = 0usize;
    let mut niter = 0usize;
    let mut frozen = false;

    while niter < maxiter && !frozen {
        niter += 1;

        // Build score table
        let table0 = filltable_maxsim(x, y, rcut);

        // Save old alignment
        let ali_old = ali.clone();

        // NW alignment
        let (_score, new_ali) = nw_maxsim(nx, ny, &table0);
        ali = new_ali;

        // Check frozen
        frozen = ali == ali_old;

        // Superpose
        let (u, t, new_lali, new_rms) = getut(nx, &ali, x, y);
        lali = new_lali;
        rms = new_rms;
        transrotate(x, &u, &t);
    }

    FitzResult {
        ali,
        rms,
        lali,
        niter,
    }
}
