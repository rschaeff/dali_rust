use std::fs;
use std::collections::HashMap;

use dali_core::parsi::{run_parsi, dowork_parsi, ParsiCd1Cache};
use dali_core::ParsiHit;

const FIXTURE_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"), "/../..",
    "/validation/fixtures/structures"
);
const PARSI_OUTPUT_5: &str = concat!(
    env!("CARGO_MANIFEST_DIR"), "/../..",
    "/validation/fixtures/parsi/parsi_output.txt"
);
const EXPANDED_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"), "/../..",
    "/validation/corpus_expansion/ground_truth_18/structures"
);
const PARSI_OUTPUT_18: &str = concat!(
    env!("CARGO_MANIFEST_DIR"), "/../..",
    "/validation/corpus_expansion/ground_truth_18/parsi/parsi_output.txt"
);

/// Parsed reference refine line.
#[derive(Debug, Clone)]
struct RefineRef {
    cd1cd2: String,
    idom: usize,
    score: i32,
    nseg: usize,
    ranges: Vec<i32>,
}

/// Parse refine output lines from ground truth file.
fn parse_refine_gt(path: &str) -> Vec<RefineRef> {
    let content = fs::read_to_string(path).unwrap();
    let mut entries = Vec::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || !line.contains("refine") {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            continue;
        }
        let cd1cd2 = parts[0].replace("refine", "");
        let idom: usize = parts[1].parse().unwrap_or(0);
        let score: i32 = parts[2].parse().unwrap_or(0);
        let nseg: usize = parts[3].parse().unwrap_or(0);
        let ranges: Vec<i32> = parts[4..].iter()
            .filter_map(|s| s.parse().ok())
            .collect();

        entries.push(RefineRef { cd1cd2, idom, score, nseg, ranges });
    }

    entries
}

/// Convert ParsiHit to refine-like format for comparison.
fn hit_to_refine(hit: &ParsiHit) -> RefineRef {
    let cd1cd2 = format!("{}{}", hit.cd1, hit.cd2);
    let mut ranges = Vec::new();
    for &(a1, a2) in &hit.ranges_cd1 {
        ranges.push(a1);
        ranges.push(a2);
    }
    for &(b1, b2) in &hit.ranges_cd2 {
        ranges.push(b1);
        ranges.push(b2);
    }
    RefineRef {
        cd1cd2,
        idom: hit.idom,
        score: hit.score,
        nseg: hit.ranges_cd1.len(),
        ranges,
    }
}

/// Compare results against ground truth with score tolerance.
/// Returns (matched, total_ref, mismatches).
fn validate_parsi(gt: &[RefineRef], results: &[RefineRef], score_tol: i32) -> (usize, usize, Vec<String>) {
    // Group by (cd1cd2, idom)
    let mut ref_map: HashMap<(String, usize), Vec<&RefineRef>> = HashMap::new();
    for r in gt {
        ref_map.entry((r.cd1cd2.clone(), r.idom)).or_default().push(r);
    }

    let mut cand_map: HashMap<(String, usize), Vec<&RefineRef>> = HashMap::new();
    for c in results {
        cand_map.entry((c.cd1cd2.clone(), c.idom)).or_default().push(c);
    }

    let mut matched = 0;
    let mut mismatches = Vec::new();

    for (key, ref_group) in &ref_map {
        let cand_group = cand_map.get(key).map(|v| v.as_slice()).unwrap_or(&[]);

        for r in ref_group {
            let mut found = false;
            // Exact match first
            for c in cand_group {
                if r.score == c.score && r.nseg == c.nseg && r.ranges == c.ranges {
                    found = true;
                    break;
                }
            }
            if !found {
                // Score tolerance match
                for c in cand_group {
                    if (r.score - c.score).abs() <= score_tol && r.nseg == c.nseg {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                matched += 1;
            } else {
                if mismatches.len() < 20 {
                    mismatches.push(format!(
                        "{} idom={} score={} nseg={}",
                        key.0, key.1, r.score, r.nseg
                    ));
                }
            }
        }
    }

    (matched, gt.len(), mismatches)
}

#[test]
fn test_parsi_5_structures() {
    let gt = parse_refine_gt(PARSI_OUTPUT_5);
    assert_eq!(gt.len(), 435, "Expected 435 ground truth refine lines for 5-struct corpus");

    let structures: Vec<String> = vec![
        "101mA", "1a00A", "1a87A", "1allA", "1binA",
    ].into_iter().map(String::from).collect();

    let hits = run_parsi(&structures, FIXTURE_DIR, None);
    let results: Vec<RefineRef> = hits.iter().map(hit_to_refine).collect();

    eprintln!("PARSI 5-struct: {} result lines vs {} reference lines",
              results.len(), gt.len());

    let (matched, total, mismatches) = validate_parsi(&gt, &results, 100);
    let pct = matched as f64 / total as f64 * 100.0;
    eprintln!("PARSI 5-struct: {}/{} matched ({:.1}%)", matched, total, pct);
    if !mismatches.is_empty() {
        for m in &mismatches {
            eprintln!("  MISSING: {}", m);
        }
    }

    // Expected: ~87% (379/435) due to search divergence
    assert!(
        matched >= 300,
        "PARSI 5-struct: only {}/{} matched ({:.1}%), expected >=300 (~69%)",
        matched, total, pct
    );
}

#[test]
fn test_parsi_18_structures() {
    let gt = parse_refine_gt(PARSI_OUTPUT_18);
    assert_eq!(gt.len(), 13515, "Expected 13515 ground truth refine lines for 18-struct corpus");

    let structures: Vec<String> = vec![
        "1a7sA", "1aiwA", "1a25A", "1a12A", "1f3uA", "1a04A",
        "1bbhA", "1b3qA", "1a17A", "1miwA", "1aopA", "1a8lA",
        "1a6qA", "1a06A", "1b3oB", "1a0cA", "1a4iA", "1bcoA",
    ].into_iter().map(String::from).collect();

    let hits = run_parsi(&structures, EXPANDED_DIR, None);
    let results: Vec<RefineRef> = hits.iter().map(hit_to_refine).collect();

    eprintln!("PARSI 18-struct: {} result lines vs {} reference lines",
              results.len(), gt.len());

    let (matched, total, mismatches) = validate_parsi(&gt, &results, 100);
    let pct = matched as f64 / total as f64 * 100.0;
    eprintln!("PARSI 18-struct: {}/{} matched ({:.1}%)", matched, total, pct);
    if !mismatches.is_empty() {
        for m in &mismatches[..mismatches.len().min(10)] {
            eprintln!("  MISSING: {}", m);
        }
    }

    // Expected: ~78% (10533/13515) due to search divergence
    assert!(
        matched >= 8000,
        "PARSI 18-struct: only {}/{} matched ({:.1}%), expected >=8000 (~59%)",
        matched, total, pct
    );
}
