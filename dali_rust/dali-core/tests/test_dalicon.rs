use std::fs;

use dali_core::dalicon::{parse_dalicon_input, run_dalicon};
use dali_core::store::ProteinStore;

/// Path to fixture directories.
const FIXTURE_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"), "/../..", "/validation/fixtures/structures"
);
const DALICON_INPUT_5: &str = concat!(
    env!("CARGO_MANIFEST_DIR"), "/../..",
    "/validation/fixtures/dalicon/wolf/dalicon_input.txt"
);
const DALICON_OUTPUT_5: &str = concat!(
    env!("CARGO_MANIFEST_DIR"), "/../..",
    "/validation/fixtures/dalicon/wolf/dalicon_output.txt"
);
const EXPANDED_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"), "/../..",
    "/validation/corpus_expansion/ground_truth_18/structures"
);
const DALICON_INPUT_18: &str = concat!(
    env!("CARGO_MANIFEST_DIR"), "/../..",
    "/validation/corpus_expansion/ground_truth_18/dalicon/wolf/dalicon_input_from_dp.txt"
);
const DALICON_OUTPUT_18: &str = concat!(
    env!("CARGO_MANIFEST_DIR"), "/../..",
    "/validation/corpus_expansion/ground_truth_18/dalicon/wolf/dalicon_output.txt"
);

// 18-structure groups (6 each)
const GROUP_A: &[&str] = &["1a7sA", "1aiwA", "1a25A", "1a12A", "1f3uA", "1a04A"];
const GROUP_B: &[&str] = &["1bbhA", "1b3qA", "1a17A", "1miwA", "1aopA", "1a8lA"];
const GROUP_C: &[&str] = &["1a6qA", "1a06A", "1b3oB", "1a0cA", "1a4iA", "1bcoA"];

/// Parsed WOLFITZ ground truth entry.
#[derive(Debug)]
struct WolfitzRef {
    cd1: String,
    cd2: String,
    nblock: usize,
    raw_values: Vec<i64>,
}

/// Parse WOLFITZ ground truth file.
fn parse_wolfitz_gt(path: &str) -> Vec<WolfitzRef> {
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
        let cd1 = pair[..5].to_string();
        let cd2 = pair[5..].to_string();
        let nblock: usize = parts[2].parse().unwrap();
        let raw_values: Vec<i64> = parts[3..]
            .iter()
            .map(|s| s.parse().unwrap())
            .collect();

        entries.push(WolfitzRef {
            cd1,
            cd2,
            nblock,
            raw_values,
        });
    }

    entries
}

/// Convert DALICON result to raw values in WOLFITZ format.
fn result_to_raw_values(blocks: &[dali_core::AlignmentBlock]) -> Vec<i64> {
    let mut vals = Vec::new();
    for b in blocks {
        vals.push(b.l1 as i64);
        vals.push(b.r1 as i64);
    }
    for b in blocks {
        vals.push(b.l2 as i64);
        vals.push(b.r2 as i64);
    }
    vals
}

/// Validate DALICON results against ground truth.
fn validate_dalicon(
    gt: &[WolfitzRef],
    results: &[(String, String, usize, Vec<dali_core::AlignmentBlock>)],
) -> (usize, usize, Vec<String>) {
    let mut matched = 0;
    let mut mismatched = Vec::new();

    let result_map: std::collections::HashMap<(&str, &str), &(String, String, usize, Vec<dali_core::AlignmentBlock>)> =
        results.iter().map(|r| ((r.0.as_str(), r.1.as_str()), r)).collect();

    for entry in gt {
        let key = (entry.cd1.as_str(), entry.cd2.as_str());
        match result_map.get(&key) {
            Some(result) => {
                let mut issues = Vec::new();
                let result_nblock = result.2;

                if entry.nblock != result_nblock {
                    issues.push(format!(
                        "nblock: ref={} got={}",
                        entry.nblock, result_nblock
                    ));
                }

                let result_raw = result_to_raw_values(&result.3);
                if entry.raw_values != result_raw {
                    for (idx, (rv, cv)) in entry.raw_values.iter()
                        .zip(result_raw.iter())
                        .enumerate()
                    {
                        if rv != cv {
                            issues.push(format!("val[{}]: ref={} got={}", idx, rv, cv));
                            break;
                        }
                    }
                    if entry.raw_values.len() != result_raw.len() {
                        issues.push(format!(
                            "val count: ref={} got={}",
                            entry.raw_values.len(),
                            result_raw.len()
                        ));
                    }
                }

                if issues.is_empty() {
                    matched += 1;
                } else {
                    mismatched.push(format!(
                        "{}{}: {}",
                        entry.cd1, entry.cd2, issues.join("; ")
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

#[test]
fn test_dalicon_5_structures() {
    let gt = parse_wolfitz_gt(DALICON_OUTPUT_5);
    assert_eq!(gt.len(), 10, "Expected 10 ground truth entries for 5-struct corpus");

    let input_content = fs::read_to_string(DALICON_INPUT_5).unwrap();
    let records = parse_dalicon_input(&input_content);

    let store = ProteinStore::new(FIXTURE_DIR);
    let results = run_dalicon(&records, &store);

    let (matched, total, mismatched) = validate_dalicon(&gt, &results);

    eprintln!("DALICON 5-struct: {}/{} matched", matched, total);
    if !mismatched.is_empty() {
        for m in &mismatched {
            eprintln!("  {}", m);
        }
    }
    assert_eq!(
        matched, 10,
        "DALICON 5-struct: {}/10 (mismatches: {:?})",
        matched, mismatched
    );
}

fn run_dalicon_18_group(codes: &[&str], label: &str) {
    let all_gt = parse_wolfitz_gt(DALICON_OUTPUT_18);
    let group_gt: Vec<&WolfitzRef> = all_gt.iter()
        .filter(|e| codes.contains(&e.cd1.as_str()) && codes.contains(&e.cd2.as_str()))
        .collect();
    let gt_count = group_gt.len();

    if gt_count == 0 {
        eprintln!("DALICON {}: no ground truth entries for this group", label);
        return;
    }

    let input_content = fs::read_to_string(DALICON_INPUT_18).unwrap();
    let all_records = parse_dalicon_input(&input_content);

    // Filter records to only those involving group structures
    let records: Vec<_> = all_records.into_iter()
        .filter(|r| codes.contains(&r.cd1.as_str()) && codes.contains(&r.cd2.as_str()))
        .collect();

    let store = ProteinStore::new(EXPANDED_DIR);
    let results = run_dalicon(&records, &store);

    // Validate against filtered ground truth
    let mut matched = 0;
    let mut mismatched = Vec::new();
    let result_map: std::collections::HashMap<(&str, &str), &(String, String, usize, Vec<dali_core::AlignmentBlock>)> =
        results.iter().map(|r| ((r.0.as_str(), r.1.as_str()), r)).collect();

    for entry in &group_gt {
        let key = (entry.cd1.as_str(), entry.cd2.as_str());
        match result_map.get(&key) {
            Some(result) => {
                let result_raw = result_to_raw_values(&result.3);
                if entry.nblock == result.2 && entry.raw_values == result_raw {
                    matched += 1;
                } else {
                    mismatched.push(format!(
                        "{}{}: nblock ref={} got={}",
                        entry.cd1, entry.cd2, entry.nblock, result.2
                    ));
                }
            }
            None => {
                mismatched.push(format!("{}{}: MISSING", entry.cd1, entry.cd2));
            }
        }
    }

    eprintln!("DALICON {}: {}/{} matched", label, matched, gt_count);
    if !mismatched.is_empty() {
        for m in &mismatched {
            eprintln!("  {}", m);
        }
    }
    assert_eq!(
        matched, gt_count,
        "DALICON {}: {}/{} (mismatches: {:?})",
        label, matched, gt_count, mismatched
    );
}

#[test]
fn test_dalicon_18_group_a() {
    run_dalicon_18_group(GROUP_A, "18-group-A");
}

#[test]
fn test_dalicon_18_group_b() {
    run_dalicon_18_group(GROUP_B, "18-group-B");
}

#[test]
fn test_dalicon_18_group_c() {
    run_dalicon_18_group(GROUP_C, "18-group-C");
}
