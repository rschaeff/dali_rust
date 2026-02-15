use std::fs;
use std::collections::HashMap;

use dali_core::parsi::run_parsi;
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

// 18-structure groups (6 each)
const GROUP_A: &[&str] = &["1a7sA", "1aiwA", "1a25A", "1a12A", "1f3uA", "1a04A"];
const GROUP_B: &[&str] = &["1bbhA", "1b3qA", "1a17A", "1miwA", "1aopA", "1a8lA"];
const GROUP_C: &[&str] = &["1a6qA", "1a06A", "1b3oB", "1a0cA", "1a4iA", "1bcoA"];

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
fn validate_parsi(gt: &[RefineRef], results: &[RefineRef], score_tol: i32) -> (usize, usize, Vec<String>) {
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
            for c in cand_group {
                if r.score == c.score && r.nseg == c.nseg && r.ranges == c.ranges {
                    found = true;
                    break;
                }
            }
            if !found {
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

/// Filter ground truth to entries where both cd1 and cd2 are in codes.
fn filter_refine_gt(gt: &[RefineRef], codes: &[&str]) -> Vec<RefineRef> {
    gt.iter()
        .filter(|r| {
            if r.cd1cd2.len() != 10 { return false; }
            let cd1 = &r.cd1cd2[..5];
            let cd2 = &r.cd1cd2[5..];
            codes.contains(&cd1) && codes.contains(&cd2)
        })
        .cloned()
        .collect()
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

    assert!(
        matched >= 300,
        "PARSI 5-struct: only {}/{} matched ({:.1}%), expected >=300 (~69%)",
        matched, total, pct
    );
}

fn run_parsi_18_group(codes: &[&str], label: &str) {
    let all_gt = parse_refine_gt(PARSI_OUTPUT_18);
    let group_gt = filter_refine_gt(&all_gt, codes);
    let gt_count = group_gt.len();

    if gt_count == 0 {
        eprintln!("PARSI {}: no ground truth entries for this group", label);
        return;
    }

    let structures: Vec<String> = codes.iter().map(|s| s.to_string()).collect();
    let hits = run_parsi(&structures, EXPANDED_DIR, None);
    let results: Vec<RefineRef> = hits.iter().map(hit_to_refine).collect();

    eprintln!("PARSI {}: {} result lines vs {} reference lines",
              label, results.len(), gt_count);

    let (matched, total, mismatches) = validate_parsi(&group_gt, &results, 100);
    let pct = matched as f64 / total.max(1) as f64 * 100.0;
    eprintln!("PARSI {}: {}/{} matched ({:.1}%)", label, matched, total, pct);
    if !mismatches.is_empty() {
        for m in &mismatches[..mismatches.len().min(10)] {
            eprintln!("  MISSING: {}", m);
        }
    }

    // Expect at least 59% match (search divergence)
    assert!(
        pct >= 59.0 || total < 10,
        "PARSI {}: only {}/{} matched ({:.1}%), expected >=59%",
        label, matched, total, pct
    );
}

#[test]
fn test_parsi_18_group_a() {
    run_parsi_18_group(GROUP_A, "18-group-A");
}

#[test]
fn test_parsi_18_group_b() {
    run_parsi_18_group(GROUP_B, "18-group-B");
}

#[test]
fn test_parsi_18_group_c() {
    run_parsi_18_group(GROUP_C, "18-group-C");
}
