//! search_database / iterative_search in both directions, as dali.pl pairwise does.

use dali_core::pipeline::{iterative_search, search_database, swap_orientation};
use dali_core::{AlignmentBlock, DccpEntry, ProteinStore};

const DAT_DIR_5: &str = concat!(
    env!("CARGO_MANIFEST_DIR"), "/../..", "/validation/fixtures/structures"
);

const PAIRS: &[(&str, &str)] = &[
    ("101mA", "1a00A"),
    ("1a87A", "1allA"),
    ("1a87A", "101mA"),
    ("1allA", "1binA"),
];

fn best(v: &[DccpEntry]) -> f64 {
    v.iter().map(|e| e.zscore).fold(f64::MIN, f64::max)
}

#[test]
fn test_swap_orientation() {
    let e = DccpEntry {
        cd1: "t".into(), cd2: "q".into(), score: 1.0, zscore: 5.0, rmsd: 2.0,
        blocks: vec![
            AlignmentBlock { l1: 1, r1: 5, l2: 40, r2: 44 },
            AlignmentBlock { l1: 10, r1: 12, l2: 3, r2: 5 },
        ],
    };
    let s = swap_orientation(&e);
    assert_eq!((s.cd1.as_str(), s.cd2.as_str()), ("q", "t"));
    assert_eq!(s.blocks, vec![
        AlignmentBlock { l1: 3, r1: 5, l2: 10, r2: 12 },
        AlignmentBlock { l1: 40, r1: 44, l2: 1, r2: 5 },
    ]);
    assert_eq!(swap_orientation(&s).blocks.len(), 2);
}

/// Every hit is reported in query orientation, with block ranges inside each protein.
#[test]
fn test_bothways_hits_are_query_oriented() {
    let store = ProteinStore::new(DAT_DIR_5);
    for &(q, t) in PAIRS {
        let nq = store.get_protein(q).unwrap().nres as u32;
        let nt = store.get_protein(t).unwrap().nres as u32;
        for e in search_database(q, &[t], &store, 2.0, false, true) {
            assert_eq!((e.cd1.as_str(), e.cd2.as_str()), (q, t));
            for b in &e.blocks {
                assert!(b.l1 >= 1 && b.r1 <= nq && b.l2 >= 1 && b.r2 <= nt,
                        "{} {}: block {:?} outside {}x{}", q, t, b, nq, nt);
            }
        }
    }
}

/// Searching both ways can only add candidates, and the result is symmetric in Z:
/// q->t bothways and t->q bothways see the same four paths.
#[test]
fn test_bothways_superset_and_symmetric() {
    let store = ProteinStore::new(DAT_DIR_5);
    for &(q, t) in PAIRS {
        let fwd = search_database(q, &[t], &store, 2.0, false, false);
        let both = search_database(q, &[t], &store, 2.0, false, true);
        let both_rev = search_database(t, &[q], &store, 2.0, false, true);
        assert!(best(&both) >= best(&fwd), "{} {}: {} < {}", q, t, best(&both), best(&fwd));
        assert!((best(&both) - best(&both_rev)).abs() < 1e-9,
                "{} {}: {} vs {}", q, t, best(&both), best(&both_rev));
    }
}

#[test]
fn test_iterative_search_bothways_round0_matches_search() {
    // iterative_search writes masked queries into the store, so work on a copy.
    let dir = std::env::temp_dir().join(format!("dali_iter_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    for c in ["1a87A", "1allA"] {
        std::fs::copy(format!("{}/{}.dat", DAT_DIR_5, c), dir.join(format!("{}.dat", c))).unwrap();
    }
    let store = ProteinStore::new(&dir);
    let hits = iterative_search("1a87A", &["1allA"], &store, 20, 2.0, 5, 3, false, true);
    let both = search_database("1a87A", &["1allA"], &store, 2.0, false, true);
    assert!(!hits.is_empty());
    assert!((hits[0].zscore - best(&both)).abs() < 1e-9);
    let _ = std::fs::remove_dir_all(&dir);
}
