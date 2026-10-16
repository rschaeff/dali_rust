
/// DP distance-difference threshold.
const D0: f32 = 0.20;

// The DP functions below are in REAL (f32), as in comparemodules.f; the
// printed DCCP score (f8.1) is sensitive to double-precision accumulation.

/// exp(-0.0025*i*i), i = 0..100, exactly as DaliLite computes it.
///
/// gfortran -O3 vectorizes the weight loops in dpweights (comparemodules.f)
/// and gagaweights (gagatool.f) into calls to glibc libmvec's
/// `_ZGVbN4v_expf`, which differs from scalar expf in 58 of the 101 entries
/// (up to ~25 ulp in the tail). DP and DALICON scores depend on these bits,
/// so the table is the one the reference binary produces (glibc 2.35 libmvec
/// on x86_64; the DP and DALICON loops give identical tables).
pub const ENVELOPE_WEIGHT_BITS: [u32; 101] = [
    0x3F800000, 0x3F7F5C5D, 0x3F7D73E8, 0x3F7A4DE8, 0x3F75F64C, 0x3F707D60,
    0x3F69F764, 0x3F627C15, 0x3F5A2618, 0x3F51125E, 0x3F475F7E, 0x3F3D2D0C,
    0x3F329AEA, 0x3F27C8B6, 0x3F1CD516, 0x3F11DD4A, 0x3F06FCA3, 0x3EF8984A,
    0x3EE3C470, 0x3ECFA4D7, 0x3EBC5AB4, 0x3EAA0126, 0x3E98AD50, 0x3E886E88,
    0x3E729D32, 0x3E56A45E, 0x3E3CF2B2, 0x3E258009, 0x3E103D30, 0x3DFA2A22,
    0x3DD7DB8F, 0x3DB95396, 0x3D9E51F4, 0x3D86931D, 0x3D63A3E0, 0x3D3F928F,
    0x3D206A55, 0x3D05A7D3, 0x3CDD9BD6, 0x3CB6CDC7, 0x3C960AB0, 0x3C7512F9,
    0x3C472686, 0x3C210671, 0x3C018C95, 0x3BCF6980, 0x3BA5355B, 0x3B82EF95,
    0x3B4E82E4, 0x3B220ABF, 0x3AFD074B, 0x3AC49117, 0x3A97F149, 0x3A69B9DE,
    0x3A32DE48, 0x3A083411, 0x39CE655B, 0x399B99EF, 0x396971C1, 0x392E3DEA,
    0x39016799, 0x38BF4079, 0x388C9FCC, 0x384DC407, 0x3815CA63, 0x37D8FF52,
    0x379C6528, 0x37604F9C, 0x37200EAC, 0x36E3477B, 0x36A08FD9, 0x3661BA22,
    0x361DE0F9, 0x35DBBF1A, 0x35982A8C, 0x3551B025, 0x350FC1AF, 0x34C420D8,
    0x34851F57, 0x3433D014, 0x33F1AAED, 0x33A1970D, 0x335703FA, 0x330E56B9,
    0x32BB83C8, 0x3275CBEA, 0x32204B0D, 0x31D005FF, 0x31864F5D, 0x312C91BD,
    0x30DC9EF0, 0x308C52AA, 0x30319C4A, 0x2FDFAF75, 0x2F8C275A, 0x2F2EC162,
    0x2ED8D022, 0x2E85D348, 0x2E246163, 0x2DC8E794, 0x2D7451BD,
];

/// The DaliLite envelope weight table as f32.
pub fn envelope_weights() -> [f32; 101] {
    ENVELOPE_WEIGHT_BITS.map(f32::from_bits)
}

/// Compute Gaussian envelope weights for DP scoring.
///
/// wght[i] = exp(-(i/radius)^2) for i=0..100, radius 20: see
/// [`ENVELOPE_WEIGHT_BITS`].
pub fn dpweights() -> [f32; 101] {
    envelope_weights()
}

/// Distance-based scoring function for DP.
///
/// Translated from comparemodules.f dpscorefun(); a-b and a+b are integer*2.
///
/// # Arguments
/// * `a`, `b` - int16 distances (dist * 10)
/// * `wght` - weight array from dpweights()
pub fn dpscorefun(a: i16, b: i16, wght: &[f32; 101]) -> f32 {
    let x = a.wrapping_sub(b).wrapping_abs() as f32 / 10.0f32;
    let y = a.wrapping_add(b) as f32 / 20.0f32;

    if y > 100.0f32 {
        return 0.0;
    }

    let iy = super::rounding::nint_f32(y).clamp(0, 100) as usize;

    if y > 0.0f32 {
        wght[iy] * (D0 - x / y)
    } else {
        wght[iy] * D0
    }
}

/// Compute Z-score from domain lengths and raw score.
///
/// Uses cubic polynomial calibration.
/// Translated from comparemodules.f zscore().
pub fn zscore_func(l1: usize, l2: usize, score: f32) -> f32 {
    let n12 = ((l1 * l2) as i32 as f32).sqrt();
    let x = n12.min(400.0f32);
    let mut mean = 7.9494f32 + 0.70852f32 * x + 2.5895e-4f32 * x * x - 1.9156e-6f32 * x * x * x;
    if n12 > 400.0f32 {
        mean += (n12 - 400.0f32) * 1.0f32;
    }
    let sigma = 0.50f32 * mean;
    (score - mean) / sigma.max(1.0f32)
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
    wght: &[f32; 101],
) -> f32 {
    // Collect aligned residue indices (0-based)
    let aligned: Vec<usize> = (0..nres1).filter(|&i| ali1[i] != 0).collect();

    let mut tots = 0.0f32;
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
        assert!((w[0] - 1.0).abs() < 1e-7);
        assert!(w[100] > 0.0);
        assert!(w[100] < w[0]);
        // w[20] = exp(-1) = 0.367879...
        assert!((w[20] - (-1.0f32).exp()).abs() < 1e-6);
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
        assert!((w[0] - 1.0).abs() < 1e-7);
        // w[10] = exp(-(10/10)^2) = exp(-1)
        assert!((w[10] - (-1.0f64).exp()).abs() < 1e-10);
    }
}
