use super::rounding::nint;

/// DP envelope radius constant.
const ENVELOPE_RADIUS: f64 = 20.0;

/// DP distance-difference threshold.
const D0: f64 = 0.20;

/// Compute Gaussian envelope weights for DP scoring.
///
/// wght[i] = exp(-(i/radius)^2) for i=0..100
pub fn dpweights() -> [f64; 101] {
    let mut wght = [0.0f64; 101];
    let x = 1.0 / (ENVELOPE_RADIUS * ENVELOPE_RADIUS);
    for i in 0..101 {
        let fi = i as f64;
        wght[i] = (-x * fi * fi).exp();
    }
    wght
}

/// Distance-based scoring function for DP.
///
/// Translated from comparemodules.f dpscorefun().
///
/// # Arguments
/// * `a`, `b` - int16 distances (dist * 10)
/// * `wght` - weight array from dpweights()
pub fn dpscorefun(a: i16, b: i16, wght: &[f64; 101]) -> f64 {
    let ai = a as i32;
    let bi = b as i32;
    let x = (ai - bi).unsigned_abs() as f64 / 10.0;
    let y = (ai + bi) as f64 / 20.0;

    if y > 100.0 {
        return 0.0;
    }

    let mut iy = nint(y);
    if iy < 0 {
        iy = 0;
    }
    if iy > 100 {
        iy = 100;
    }

    if y > 0.0 {
        wght[iy as usize] * (D0 - x / y)
    } else {
        wght[iy as usize] * D0
    }
}

/// Compute Z-score from domain lengths and raw score.
///
/// Uses cubic polynomial calibration.
/// Translated from comparemodules.f zscore().
pub fn zscore_func(l1: usize, l2: usize, score: f64) -> f64 {
    let n12 = ((l1 * l2) as f64).sqrt();
    let x = n12.min(400.0);
    let mut mean = 7.9494 + 0.70852 * x + 2.5895e-4 * x * x - 1.9156e-6 * x * x * x;
    if n12 > 400.0 {
        mean += (n12 - 400.0) * 1.0;
    }
    let sigma = 0.50 * mean;
    (score - mean) / sigma.max(1.0)
}

/// Compute total alignment score.
///
/// Double loop over all pairs of aligned residues.
/// Translated from comparemodules.f totscore().
///
/// # Arguments
/// * `ali1` - alignment array (length nres1). ali1[i]=j (1-based target) means
///   residue i aligned to residue j. 0 means unaligned.
/// * `nres1` - number of residues in protein 1
/// * `d1` - (nres1, nres1) distance matrix for protein 1
/// * `d2` - (nres2, nres2) distance matrix for protein 2
/// * `wght` - weight array
pub fn totscore(
    ali1: &[i32],
    nres1: usize,
    d1: &ndarray::Array2<i16>,
    d2: &ndarray::Array2<i16>,
    wght: &[f64; 101],
) -> f64 {
    // Collect aligned residue indices (0-based)
    let aligned: Vec<usize> = (0..nres1).filter(|&i| ali1[i] != 0).collect();

    let mut tots = 0.0;
    for &k in &aligned {
        for &l in &aligned {
            let q = (ali1[k].unsigned_abs() - 1) as usize; // 1-based target to 0-based
            let r = (ali1[l].unsigned_abs() - 1) as usize;
            tots += dpscorefun(d1[[k, l]], d2[[q, r]], wght);
        }
    }

    tots
}

/// Compute pairwise CA distance matrix as int16 (dist * 10).
///
/// Uses f32 intermediate arithmetic to match Fortran's single-precision.
/// Translated from comparemodules.f dpgetdist().
pub fn dpgetdist(ca: &ndarray::Array2<f64>, nres: usize) -> ndarray::Array2<i16> {
    let mut d = ndarray::Array2::<i16>::zeros((nres, nres));
    let ten: f32 = 10.0;

    for i in 0..nres {
        let xi = ca[[0, i]] as f32;
        let yi = ca[[1, i]] as f32;
        let zi = ca[[2, i]] as f32;
        for j in 0..i {
            let dx = xi - ca[[0, j]] as f32;
            let dy = yi - ca[[1, j]] as f32;
            let dz = zi - ca[[2, j]] as f32;
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            let x = super::rounding::nint_f32(ten * dist);
            let x = x.min(32767);
            d[[i, j]] = x as i16;
            d[[j, i]] = x as i16;
        }
    }

    d
}

/// Gaussian weights for DALICON scoring.
///
/// gagaweights(i) = exp(-(i/10)^2) for i=0..100
pub fn gagaweights() -> [f64; 101] {
    let mut wght = [0.0f64; 101];
    let x = 1.0 / 100.0; // (1/10)^2
    for i in 0..101 {
        let fi = i as f64;
        wght[i] = (-x * fi * fi).exp();
    }
    wght
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dpweights() {
        let w = dpweights();
        assert!((w[0] - 1.0).abs() < 1e-15);
        assert!(w[100] > 0.0);
        assert!(w[100] < w[0]);
        // w[20] = exp(-1) = 0.367879...
        assert!((w[20] - (-1.0f64).exp()).abs() < 1e-10);
    }

    #[test]
    fn test_dpscorefun_identical() {
        let w = dpweights();
        // Identical distances: x=0, score = w[iy] * D0
        let s = dpscorefun(50, 50, &w);
        // y = 100/20 = 5.0, iy=5, w[5]=exp(-25/400)=exp(-0.0625)
        // s = w[5] * (0.20 - 0/5) = w[5] * 0.20
        assert!(s > 0.0);
    }

    #[test]
    fn test_zscore_func() {
        // Very similar structures should have high z-score
        let z = zscore_func(100, 100, 500.0);
        assert!(z > 2.0);

        // Random score should be near zero
        let z = zscore_func(100, 100, 80.0);
        assert!(z.abs() < 2.0);
    }

    #[test]
    fn test_gagaweights() {
        let w = gagaweights();
        assert!((w[0] - 1.0).abs() < 1e-15);
        // w[10] = exp(-(10/10)^2) = exp(-1)
        assert!((w[10] - (-1.0f64).exp()).abs() < 1e-10);
    }
}
