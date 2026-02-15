use dali_core::numerics::scoring::{dpweights, dpscorefun, zscore_func};

#[test]
fn test_dpscorefun_values() {
    let w = dpweights();

    // Reference values from Python
    let cases: Vec<(i16, i16, f64)> = vec![
        (38, 38, 0.1921578878),
        (50, 100, -0.3976671015),
        (0, 0, 0.2000000000),
        (200, 200, 0.0735758882),
        (100, 150, -0.1310812509),
    ];

    for (a, b, expected) in cases {
        let result = dpscorefun(a, b, &w);
        assert!(
            (result - expected).abs() < 1e-8,
            "dpscorefun({}, {}) = {}, expected {}",
            a, b, result, expected
        );
    }
}

#[test]
fn test_zscore_values() {
    let cases: Vec<(usize, usize, f64, f64)> = vec![
        (100, 100, 500.0, 10.5825256400),
        (154, 141, 200.0, 1.5762800414),
        (50, 50, 50.0, 0.2839745497),
    ];

    for (l1, l2, score, expected) in cases {
        let result = zscore_func(l1, l2, score);
        assert!(
            (result - expected).abs() < 1e-8,
            "zscore({}, {}, {}) = {}, expected {}",
            l1, l2, score, result, expected
        );
    }
}
