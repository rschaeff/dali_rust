use std::fs;

use dali_core::types::AlignmentBlock;
use dali_core::dp::{run_dp, ZCUT_DEFAULT};
use dali_core::store::ProteinStore;

/// Path to fixture directories.
const FIXTURE_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"), "/../..", "/validation/fixtures/structures"
);
const DP_INPUT_5: &str = concat!(
    env!("CARGO_MANIFEST_DIR"), "/../..", "/validation/fixtures/dp/wolf/dp_input.txt"
);
const DP_OUTPUT_5: &str = concat!(
    env!("CARGO_MANIFEST_DIR"), "/../..", "/validation/fixtures/dp/wolf/dp_output.txt"
);
const EXPANDED_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"), "/../..", "/validation/corpus_expansion/ground_truth_18/structures"
);
const DP_INPUT_18: &str = concat!(
    env!("CARGO_MANIFEST_DIR"), "/../..",
    "/validation/corpus_expansion/ground_truth_18/wolf/wolf_output.txt"
);
const DP_OUTPUT_18: &str = concat!(
    env!("CARGO_MANIFEST_DIR"), "/../..",
    "/validation/corpus_expansion/ground_truth_18/dp/wolf/dp_output.txt"
);

/// Parsed DCCP ground truth entry.
#[derive(Debug)]
struct DccpRef {
    cd1: String,
    cd2: String,
    score: f64,
    nblock: usize,
    ranges_cd1: Vec<(i32, i32)>,
    ranges_cd2: Vec<(i32, i32)>,
}

/// Parse DCCP ground truth file.
fn parse_dccp(path: &str) -> Vec<DccpRef> {
    let content = fs::read_to_string(path).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    let mut entries = Vec::new();
    let mut i = 0;

    while i < lines.len() {
        let line = lines[i].trim();
        if !line.starts_with("DCCP") {
            i += 1;
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        let score: f64 = parts[2].parse().unwrap();
        let nblock: usize = parts[7].parse().unwrap();
        let cd1 = parts[8].to_string();
        let cd2 = parts[9].to_string();

        // Skip "alignment" line
        i += 1;
        if i < lines.len() && lines[i].contains("alignment") {
            i += 1;
        }

        // Read range values (nblock * 4 integers total)
        let mut range_values: Vec<i32> = Vec::new();
        while i < lines.len() && range_values.len() < nblock * 4 {
            let line = lines[i].trim();
            if line.starts_with("DCCP") {
                break;
            }
            for val in line.split_whitespace() {
                if let Ok(v) = val.parse::<i32>() {
                    range_values.push(v);
                }
            }
            i += 1;
            if range_values.len() >= nblock * 4 {
                break;
            }
        }

        let mut ranges_cd1 = Vec::new();
        let mut ranges_cd2 = Vec::new();
        if range_values.len() >= nblock * 4 {
            for b in 0..nblock {
                ranges_cd1.push((range_values[b * 2], range_values[b * 2 + 1]));
            }
            let offset = nblock * 2;
            for b in 0..nblock {
                ranges_cd2.push((range_values[offset + b * 2], range_values[offset + b * 2 + 1]));
            }
        }

        entries.push(DccpRef {
            cd1, cd2, score, nblock, ranges_cd1, ranges_cd2,
        });
    }

    entries
}

/// Parse WOLFITZ file into alignment tuples for run_dp().
fn parse_wolfitz_as_alignments(path: &str) -> Vec<(String, String, Vec<AlignmentBlock>)> {
    let content = fs::read_to_string(path).unwrap();
    let mut alignments = Vec::new();

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

        let raw_values: Vec<i32> = parts[3..].iter()
            .map(|s| s.parse().unwrap())
            .collect();

        let mut blocks = Vec::new();
        for b in 0..nblock {
            let l1 = raw_values[b * 2] as u32;
            let r1 = raw_values[b * 2 + 1] as u32;
            let offset = nblock * 2;
            let l2 = raw_values[offset + b * 2] as u32;
            let r2 = raw_values[offset + b * 2 + 1] as u32;
            blocks.push(AlignmentBlock { l1, r1, l2, r2 });
        }

        alignments.push((cd1, cd2, blocks));
    }

    alignments
}

/// Compare DP results against ground truth.
/// Returns (matched, total, mismatched_details).
fn validate_dp(
    gt: &[DccpRef],
    results: &[dali_core::DccpEntry],
    score_tol: f64,
) -> (usize, usize, Vec<String>) {
    let mut matched = 0;
    let mut mismatched = Vec::new();

    // Build result map keyed by (cd1, cd2)
    let result_map: std::collections::HashMap<(&str, &str), &dali_core::DccpEntry> = results
        .iter()
        .map(|r| ((r.cd1.as_str(), r.cd2.as_str()), r))
        .collect();

    for entry in gt {
        let key = (entry.cd1.as_str(), entry.cd2.as_str());
        match result_map.get(&key) {
            Some(result) => {
                let mut issues = Vec::new();

                if (entry.score - result.score).abs() > score_tol {
                    issues.push(format!(
                        "score: ref={:.1} got={:.1}",
                        entry.score, result.score
                    ));
                }

                let result_nblock = result.blocks.len();
                if entry.nblock != result_nblock {
                    issues.push(format!(
                        "nblock: ref={} got={}",
                        entry.nblock, result_nblock
                    ));
                }

                // Compare ranges
                let result_ranges_cd1: Vec<(i32, i32)> = result.blocks.iter()
                    .map(|b| (b.l1 as i32, b.r1 as i32))
                    .collect();
                let result_ranges_cd2: Vec<(i32, i32)> = result.blocks.iter()
                    .map(|b| (b.l2 as i32, b.r2 as i32))
                    .collect();

                if entry.ranges_cd1 != result_ranges_cd1 {
                    issues.push("cd1 ranges differ".to_string());
                }
                if entry.ranges_cd2 != result_ranges_cd2 {
                    issues.push("cd2 ranges differ".to_string());
                }

                if issues.is_empty() {
                    matched += 1;
                } else {
                    mismatched.push(format!("{}{}: {}", entry.cd1, entry.cd2, issues.join("; ")));
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
fn test_dp_5_structures() {
    let gt = parse_dccp(DP_OUTPUT_5);
    assert_eq!(gt.len(), 10, "Expected 10 DCCP entries for 5-struct corpus");

    let alignments = parse_wolfitz_as_alignments(DP_INPUT_5);

    let store = ProteinStore::new(FIXTURE_DIR);
    let results = run_dp(&alignments, &store, ZCUT_DEFAULT);

    let (matched, total, mismatched) = validate_dp(&gt, &results, 1.0);

    eprintln!("DP 5-struct: {}/{} matched", matched, total);
    if !mismatched.is_empty() {
        for m in &mismatched {
            eprintln!("  {}", m);
        }
    }
    assert_eq!(matched, 10, "DP 5-struct: {}/10 (mismatches: {:?})", matched, mismatched);
}

#[test]
fn test_dp_18_structures() {
    let gt = parse_dccp(DP_OUTPUT_18);
    assert_eq!(gt.len(), 24, "Expected 24 DCCP entries for 18-struct corpus");

    let alignments = parse_wolfitz_as_alignments(DP_INPUT_18);

    let store = ProteinStore::new(EXPANDED_DIR);
    let results = run_dp(&alignments, &store, ZCUT_DEFAULT);

    let (matched, total, mismatched) = validate_dp(&gt, &results, 2.0);

    eprintln!("DP 18-struct: {}/{} matched", matched, total);
    if !mismatched.is_empty() {
        for m in &mismatched {
            eprintln!("  {}", m);
        }
    }
    assert_eq!(matched, 24, "DP 18-struct: {}/24 (mismatches: {:?})", matched, mismatched);
}
