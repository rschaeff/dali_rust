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
    // cd1 ranges: l1, r1 pairs
    for b in blocks {
        vals.push(b.l1 as i64);
        vals.push(b.r1 as i64);
    }
    // cd2 ranges: l2, r2 pairs
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

    // Build result map
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
                    // Find first difference
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

#[test]
fn test_dalicon_18_structures() {
    let gt = parse_wolfitz_gt(DALICON_OUTPUT_18);
    assert_eq!(gt.len(), 24, "Expected 24 ground truth entries for 18-struct corpus");

    let input_content = fs::read_to_string(DALICON_INPUT_18).unwrap();
    let records = parse_dalicon_input(&input_content);

    let store = ProteinStore::new(EXPANDED_DIR);
    let results = run_dalicon(&records, &store);

    let (matched, total, mismatched) = validate_dalicon(&gt, &results);

    eprintln!("DALICON 18-struct: {}/{} matched", matched, total);
    if !mismatched.is_empty() {
        for m in &mismatched {
            eprintln!("  {}", m);
        }
    }
    assert_eq!(
        matched, 24,
        "DALICON 18-struct: {}/24 (mismatches: {:?})",
        matched, mismatched
    );
}
