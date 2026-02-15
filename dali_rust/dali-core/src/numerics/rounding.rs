/// Fortran-compatible rounding: round half away from zero.
///
/// Fortran's nint(0.5)=1, nint(1.5)=2, nint(2.5)=3.
/// Rust's f64::round() already rounds half away from zero, matching Fortran.
#[inline]
pub fn nint(x: f64) -> i32 {
    x.round() as i32
}

/// nint for f32 — used in distance matrix computations where
/// Fortran uses single-precision arithmetic.
#[inline]
pub fn nint_f32(x: f32) -> i32 {
    x.round() as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nint_positive_half() {
        assert_eq!(nint(0.5), 1);
        assert_eq!(nint(1.5), 2);
        assert_eq!(nint(2.5), 3);
    }

    #[test]
    fn test_nint_negative_half() {
        assert_eq!(nint(-0.5), -1);
        assert_eq!(nint(-1.5), -2);
        assert_eq!(nint(-2.5), -3);
    }

    #[test]
    fn test_nint_regular() {
        assert_eq!(nint(0.3), 0);
        assert_eq!(nint(0.7), 1);
        assert_eq!(nint(-0.3), 0);
        assert_eq!(nint(-0.7), -1);
    }
}
