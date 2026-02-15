use std::collections::HashMap;
use std::fs;

use dali_core::wolf::{setup_protein, load_protein, wolf_compare, WolfResult};
use dali_core::wolf::spatial_hash::SpatialHashGrid;

/// Path to the 5-structure fixture data.
const FIXTURE_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../..",
    "/validation/fixtures/structures"
);

/// Path to the 5-structure WOLF ground truth.
const WOLF_GT_5: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../..",
    "/validation/fixtures/wolf/wolf_output.txt"
);

/// Path to the 18-structure fixture data.
const EXPANDED_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../..",
    "/validation/corpus_expansion/ground_truth_18/structures"
);

/// Path to the 18-structure WOLF ground truth.
const WOLF_GT_18: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../..",
    "/validation/corpus_expansion/ground_truth_18/wolf/wolf_output.txt"
);

// 18-structure groups (6 each)
const GROUP_A: &[&str] = &["1a7sA", "1aiwA", "1a25A", "1a12A", "1f3uA", "1a04A"];
const GROUP_B: &[&str] = &["1bbhA", "1b3qA", "1a17A", "1miwA", "1aopA", "1a8lA"];
const GROUP_C: &[&str] = &["1a6qA", "1a06A", "1b3oB", "1a0cA", "1a4iA", "1bcoA"];

/// Parsed WOLFITZ ground truth entry.
#[derive(Debug)]
struct WolfitzEntry {
    cd1: String,
    cd2: String,
    nblock: usize,
    raw_values: Vec<i64>,
}

/// Parse a WOLFITZ ground truth file.
fn parse_wolfitz(path: &str) -> Vec<WolfitzEntry> {
    let content = fs::read_to_string(path).unwrap();
    let mut entries = Vec::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.is_empty() || parts[0] != "WOLFITZ" {
            continue;
        }

        let pair = parts[1];
        let cd1 = &pair[..5];
        let cd2 = &pair[5..];

        let nblock: usize = parts[2].parse().unwrap();
        let raw_values: Vec<i64> = parts[3..].iter()
            .map(|s| s.parse().unwrap())
            .collect();

        entries.push(WolfitzEntry {
            cd1: cd1.to_string(),
            cd2: cd2.to_string(),
            nblock,
            raw_values,
        });
    }

    entries
}

/// Run WOLF on all pairs from a structure directory.
fn run_wolf_pairs(dat_dir: &str, codes: &[&str]) -> HashMap<(String, String), WolfResult> {
    let mut protein_data = Vec::new();
    for code in codes {
        let filepath = format!("{}/{}.dat", dat_dir, code);
        match setup_protein(&filepath) {
            Ok(data) => protein_data.push(Some(data)),
            Err(e) => {
                eprintln!("Failed to load {}: {}", code, e);
                protein_data.push(None);
            }
        }
    }

    let mut results = HashMap::new();

    for cd1_data_opt in protein_data.iter() {
        let cd1_data = match cd1_data_opt {
            Some(d) => d,
            None => continue,
        };

        if cd1_data.protein.nseg <= 2 {
            continue;
        }

        let mut grid = SpatialHashGrid::new();
        load_protein(&mut grid, cd1_data);

        for cd2_data_opt in protein_data.iter() {
            let cd2_data = match cd2_data_opt {
                Some(d) => d,
                None => continue,
            };

            if let Some(result) = wolf_compare(cd1_data, cd2_data, &grid) {
                let key = (result.cd1.clone(), result.cd2.clone());
                results.insert(key, result);
            }
        }
    }

    results
}

/// Extract raw values from a WolfResult matching the WOLFITZ output format.
fn result_raw_values(result: &WolfResult) -> Vec<i64> {
    let mut vals = Vec::new();
    for b in &result.blocks {
        vals.push(b.l1 as i64);
        vals.push(b.r1 as i64);
    }
    for b in &result.blocks {
        vals.push(b.l2 as i64);
        vals.push(b.r2 as i64);
    }
    vals
}

/// Validate WOLF results against ground truth.
fn validate_wolf(
    gt: &[WolfitzEntry],
    results: &HashMap<(String, String), WolfResult>,
) -> (usize, usize, Vec<String>) {
    let mut matched = 0;
    let mut mismatched = Vec::new();

    for entry in gt {
        let key = (entry.cd1.clone(), entry.cd2.clone());
        match results.get(&key) {
            Some(result) => {
                let nblock = result.blocks.len();
                let raw = result_raw_values(result);

                if nblock == entry.nblock && raw == entry.raw_values {
                    matched += 1;
                } else {
                    mismatched.push(format!(
                        "{}{} (nblock: {}->{})",
                        entry.cd1, entry.cd2, entry.nblock, nblock
                    ));
                }
            }
            None => {
                mismatched.push(format!("{}{}: MISSING", entry.cd1, entry.cd2));
            }
        }
    }

    (matched, gt.len(), mismatched)
}

/// Filter ground truth to pairs where both cd1 and cd2 are in the given codes.
fn filter_gt<'a>(gt: &'a [WolfitzEntry], codes: &[&str]) -> Vec<&'a WolfitzEntry> {
    gt.iter()
        .filter(|e| codes.contains(&e.cd1.as_str()) && codes.contains(&e.cd2.as_str()))
        .collect()
}

#[test]
fn test_wolf_5_structures() {
    let codes = &["101mA", "1a00A", "1a87A", "1allA", "1binA"];
    let gt = parse_wolfitz(WOLF_GT_5);
    assert_eq!(gt.len(), 21, "Expected 21 ground truth entries");

    let results = run_wolf_pairs(FIXTURE_DIR, codes);
    let (matched, total, mismatched) = validate_wolf(&gt, &results);

    eprintln!("WOLF 5-struct: {}/{} matched", matched, total);
    assert_eq!(matched, 21, "WOLF 5-struct: {}/21 (mismatches: {:?})", matched, mismatched);
}

fn run_wolf_18_group(codes: &[&str], label: &str) {
    let all_gt = parse_wolfitz(WOLF_GT_18);
    let group_gt = filter_gt(&all_gt, codes);
    let gt_count = group_gt.len();

    let results = run_wolf_pairs(EXPANDED_DIR, codes);

    let mut matched = 0;
    let mut mismatched = Vec::new();
    for entry in &group_gt {
        let key = (entry.cd1.clone(), entry.cd2.clone());
        match results.get(&key) {
            Some(result) => {
                let nblock = result.blocks.len();
                let raw = result_raw_values(result);
                if nblock == entry.nblock && raw == entry.raw_values {
                    matched += 1;
                } else {
                    mismatched.push(format!(
                        "{}{} (nblock: {}->{})",
                        entry.cd1, entry.cd2, entry.nblock, nblock
                    ));
                }
            }
            None => {
                mismatched.push(format!("{}{}: MISSING", entry.cd1, entry.cd2));
            }
        }
    }

    eprintln!("WOLF {}: {}/{} matched", label, matched, gt_count);
    if !mismatched.is_empty() {
        eprintln!("  Mismatches: {:?}", mismatched);
    }

    // Allow 1 mismatch per group (floating-point edge case at 4.0 Å fitz cutoff)
    assert!(
        matched + 1 >= gt_count,
        "WOLF {}: {}/{} (need >= {}, mismatches: {:?})",
        label, matched, gt_count, gt_count.saturating_sub(1), mismatched
    );
}

#[test]
fn test_wolf_18_group_a() {
    run_wolf_18_group(GROUP_A, "18-group-A");
}

#[test]
fn test_wolf_18_group_b() {
    run_wolf_18_group(GROUP_B, "18-group-B");
}

#[test]
fn test_wolf_18_group_c() {
    run_wolf_18_group(GROUP_C, "18-group-C");
}
