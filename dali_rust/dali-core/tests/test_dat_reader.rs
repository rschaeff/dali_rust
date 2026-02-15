use dali_core::io::dat::read_dat;
use dali_core::types::{SseType, NodeType};
use dali_core::numerics::scoring::dpgetdist;

/// Path to the 5-structure fixture data.
const FIXTURE_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../..",
    "/validation/fixtures/structures"
);

#[test]
fn test_read_101ma() {
    let p = read_dat(&format!("{}/101mA.dat", FIXTURE_DIR)).unwrap();
    assert_eq!(p.code, "101mA");
    assert_eq!(p.nres, 154);
    assert_eq!(p.nseg, 5);
    assert_eq!(p.na, 5);
    assert_eq!(p.nb, 0);
    assert_eq!(p.secstr, vec![SseType::Helix; 5]);

    // Check CA coordinates
    assert!((p.ca[[0, 0]] - 24.4).abs() < 0.01);
    assert!((p.ca[[1, 0]] - 9.9).abs() < 0.01);
    assert!((p.ca[[2, 0]] - (-9.9)).abs() < 0.01);
    assert!((p.ca[[0, 153]] - 22.8).abs() < 0.01);
    assert!((p.ca[[1, 153]] - (-2.8)).abs() < 0.01);
    assert!((p.ca[[2, 153]] - 7.2).abs() < 0.01);

    // Check segments
    assert_eq!(p.segments.len(), 5);
    assert_eq!(p.segments[0].start, 5);
    assert_eq!(p.segments[0].end, 36);
    assert_eq!(p.segments[4].start, 126);
    assert_eq!(p.segments[4].end, 149);

    // Check domain tree
    assert!(!p.domain_tree.is_empty());
    assert_eq!(p.domain_tree[0].node_type, NodeType::Root);
}

#[test]
fn test_read_1a00a() {
    let p = read_dat(&format!("{}/1a00A.dat", FIXTURE_DIR)).unwrap();
    assert_eq!(p.code, "1a00A");
    assert_eq!(p.nres, 141);
    assert_eq!(p.nseg, 5);
    assert_eq!(p.na, 5);
    assert_eq!(p.nb, 0);
    assert_eq!(p.secstr, vec![SseType::Helix; 5]);
}

#[test]
fn test_read_1a87a() {
    let p = read_dat(&format!("{}/1a87A.dat", FIXTURE_DIR)).unwrap();
    assert_eq!(p.code, "1a87A");
    assert_eq!(p.nres, 297);
    assert_eq!(p.nseg, 17);
    assert_eq!(p.na, 10);
    assert_eq!(p.nb, 7);

    // Mixed SSE types: 7 E then 10 H
    assert_eq!(p.secstr[0], SseType::Strand);
    assert_eq!(p.secstr[6], SseType::Strand);
    assert_eq!(p.secstr[7], SseType::Helix);
    assert_eq!(p.secstr[16], SseType::Helix);
}

#[test]
fn test_read_1alla() {
    let p = read_dat(&format!("{}/1allA.dat", FIXTURE_DIR)).unwrap();
    assert_eq!(p.code, "1allA");
    assert_eq!(p.nres, 160);
    assert_eq!(p.nseg, 9);
    assert_eq!(p.na, 9);
    assert_eq!(p.nb, 0);
}

#[test]
fn test_read_1bina() {
    let p = read_dat(&format!("{}/1binA.dat", FIXTURE_DIR)).unwrap();
    assert_eq!(p.code, "1binA");
    assert_eq!(p.nres, 143);
    assert_eq!(p.nseg, 6);
    assert_eq!(p.na, 6);
    assert_eq!(p.nb, 0);
}

#[test]
fn test_domain_filtering_101ma() {
    let p = read_dat(&format!("{}/101mA.dat", FIXTURE_DIR)).unwrap();
    let doms = p.filtered_domains();
    assert_eq!(doms.len(), 1);
    assert_eq!(doms[0].segments, vec![(1, 154)]);
}

#[test]
fn test_domain_filtering_1a87a() {
    let p = read_dat(&format!("{}/1a87A.dat", FIXTURE_DIR)).unwrap();
    let doms = p.filtered_domains();
    // Python says: ndom=6 with segments [(1,297)], [(1,101)], [(102,297)], etc.
    assert_eq!(doms.len(), 6);
    assert_eq!(doms[0].segments, vec![(1, 297)]);
    assert_eq!(doms[1].segments, vec![(1, 101)]);
    assert_eq!(doms[2].segments, vec![(102, 297)]);
}

#[test]
fn test_domain_filtering_1alla() {
    let p = read_dat(&format!("{}/1allA.dat", FIXTURE_DIR)).unwrap();
    let doms = p.filtered_domains();
    assert_eq!(doms.len(), 2);
    assert_eq!(doms[0].segments, vec![(1, 160)]);
    assert_eq!(doms[1].segments, vec![(31, 160)]);
}

#[test]
fn test_dpgetdist_101ma() {
    let p = read_dat(&format!("{}/101mA.dat", FIXTURE_DIR)).unwrap();
    let d = dpgetdist(&p.ca, p.nres);

    // Reference values from Python
    assert_eq!(d[[0, 1]], 38);
    assert_eq!(d[[0, 10]], 158);
    assert_eq!(d[[5, 10]], 88);
    assert_eq!(d[[50, 100]], 233);
    assert_eq!(d[[153, 0]], 214);
    assert_eq!(d[[100, 50]], 233);

    // Symmetry
    assert_eq!(d[[0, 1]], d[[1, 0]]);
    assert_eq!(d[[50, 100]], d[[100, 50]]);

    // Diagonal is zero
    assert_eq!(d[[0, 0]], 0);
    assert_eq!(d[[77, 77]], 0);
}

// Test all 18 expanded corpus structures can be loaded
#[test]
fn test_load_expanded_corpus() {
    let expanded_dir = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../..",
        "/validation/corpus_expansion/ground_truth_18/structures"
    );
    let codes = [
        "1a7sA", "1aiwA", "1a25A", "1a12A", "1f3uA", "1a04A",
        "1bbhA", "1b3qA", "1a17A", "1miwA", "1aopA", "1a8lA",
        "1a6qA", "1a06A", "1b3oB", "1a0cA", "1a4iA", "1bcoA",
    ];

    for code in &codes {
        let p = read_dat(&format!("{}/{}.dat", expanded_dir, code))
            .unwrap_or_else(|e| panic!("Failed to read {}: {}", code, e));
        assert!(p.nres > 0, "{} has nres=0", code);
        assert!(p.nseg > 0, "{} has nseg=0", code);
        assert_eq!(p.ca.ncols(), p.nres, "{} CA shape mismatch", code);
        assert_eq!(p.segments.len(), p.nseg, "{} segment count mismatch", code);
        assert!(
            !p.domain_tree.is_empty(),
            "{} has no domain tree",
            code
        );
    }
}
