use ndarray::Array2;

use super::kabsch::{u3b_f32, transrotate_f32, U3bF32};
use super::nw::{filltable_maxsim, nw_maxsim};

/// Result of iterative superposition refinement.
pub struct FitzResult {
    pub ali: Vec<i32>,     // alignment, ali[i] = j (0-based) or -1
    pub rms: f64,          // final RMSD
    pub lali: usize,       // alignment length
    pub niter: usize,      // iterations performed
}

/// subfitzfast.f getut: superimpose the aligned pairs of x onto y.
///
/// Returns (U3B result, lali, rms) with rms = sqrt(ssq/max(1,lali)) in REAL.
pub fn getut_f32(ali: &[i32], x: &Array2<f32>, y: &Array2<f32>) -> (U3bF32, usize, f32) {
    let pairs: Vec<(usize, usize)> = (0..x.ncols().min(ali.len()))
        .filter(|&i| ali[i] >= 0)
        .map(|i| (i, ali[i] as usize))
        .collect();
    let lali = pairs.len();
    let mut ux = Array2::<f32>::zeros((3, lali));
    let mut uy = Array2::<f32>::zeros((3, lali));
    for (k, &(i, j)) in pairs.iter().enumerate() {
        for dim in 0..3 {
            ux[[dim, k]] = x[[dim, i]];
            uy[[dim, k]] = y[[dim, j]];
        }
    }
    let w = vec![1.0f32; lali];
    let res = u3b_f32(&w, &ux, &uy, lali, 1);
    let rms = (res.rms / lali.max(1) as f32).sqrt();
    (res, lali, rms)
}

/// Iterative superposition refinement (subfitzfast.f fitz), in REAL.
///
/// Repeatedly: score matrix -> NW alignment -> Kabsch superposition,
/// until alignment freezes or maxiter reached. The alignment starts empty
/// (Fortran zeroes ali on entry), so the first iteration never freezes on
/// a non-empty alignment.
///
/// Modifies x in-place (applies rotation/translation each iteration).
pub fn fitz(x: &mut Array2<f32>, y: &Array2<f32>, rcut: f32, maxiter: usize) -> FitzResult {
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
    let mut rms = 0.0f32;
    let mut lali = 0usize;
    let mut niter = 0usize;
    let mut frozen = false;

    while niter < maxiter && !frozen {
        niter += 1;

        let table0 = filltable_maxsim(x, y, rcut);
        let ali_old = ali.clone();
        let (_score, new_ali) = nw_maxsim(nx, ny, &table0);
        ali = new_ali;
        frozen = ali == ali_old;

        let (res, new_lali, new_rms) = getut_f32(&ali, x, y);
        lali = new_lali;
        rms = new_rms;
        transrotate_f32(x, &res.u, &res.t);
    }

    FitzResult {
        ali,
        rms: rms as f64,
        lali,
        niter,
    }
}
