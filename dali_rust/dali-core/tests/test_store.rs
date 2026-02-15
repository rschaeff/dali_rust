use dali_core::store::ProteinStore;

const FIXTURE_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../..",
    "/validation/fixtures/structures"
);

#[test]
fn test_store_load_and_cache() {
    let store = ProteinStore::new(FIXTURE_DIR);

    assert!(store.is_empty());
    assert!(!store.has_protein("101mA"));

    // Load protein
    let p = store.get_protein("101mA").unwrap();
    assert_eq!(p.nres, 154);
    assert_eq!(p.nseg, 5);

    // Should be cached now
    assert!(store.has_protein("101mA"));
    assert_eq!(store.len(), 1);

    // Load another
    let p2 = store.get_protein("1a87A").unwrap();
    assert_eq!(p2.nres, 297);
    assert_eq!(store.len(), 2);

    // Reload from cache should give same data
    let p_again = store.get_protein("101mA").unwrap();
    assert_eq!(p_again.nres, 154);
}

#[test]
fn test_store_dist_scale10() {
    let store = ProteinStore::new(FIXTURE_DIR);
    let d = store.get_dist_scale10("101mA").unwrap();

    assert_eq!(d.nres, 154);
    assert_eq!(d.data[[0, 1]], 38);
    assert_eq!(d.data[[50, 100]], 233);
}

#[test]
fn test_store_dist_scale100() {
    let store = ProteinStore::new(FIXTURE_DIR);
    let d = store.get_dist_scale100("101mA").unwrap();

    assert_eq!(d.nres, 154);
    // Scale 100 values should be ~10x the scale 10 values
    let d10 = store.get_dist_scale10("101mA").unwrap();
    let ratio = d.data[[0, 1]] as f64 / d10.data[[0, 1]] as f64;
    assert!((ratio - 10.0).abs() < 0.5, "ratio should be ~10, got {}", ratio);
}

#[test]
fn test_store_missing_file() {
    let store = ProteinStore::new(FIXTURE_DIR);
    let result = store.get_protein("nonexistent");
    assert!(result.is_err());
}
