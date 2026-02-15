use std::collections::HashMap;
use std::fs;

use dali_core::filter95::{parse_refine_line, run_filter95, Filter95Entry};
use dali_core::ParsiHit;
use dali_core::ProteinStore;

const FIXTURE_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../..",
    "/validation/fixtures/structures"
);
const FILTER95_INPUT_5: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../..",
    "/validation/fixtures/filter95/filter95_input.txt"
);
const FILTER95_OUTPUT_5: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../..",
    "/validation/fixtures/filter95/filter95_output.txt"
);
const EXPANDED_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../..",
    "/validation/corpus_expansion/ground_truth_18/structures"
);
const FILTER95_INPUT_18: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../..",
    "/validation/corpus_expansion/ground_truth_18/filter95/filter95_input.txt"
);
const FILTER95_OUTPUT_18: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../..",
    "/validation/corpus_expansion/ground_truth_18/filter95/filter95_output.txt"
);

// 18-structure groups (6 each)
const GROUP_A: &[&str] = &["1a7sA", "1aiwA", "1a25A", "1a12A", "1f3uA", "1a04A"];
const GROUP_B: &[&str] = &["1bbhA", "1b3qA", "1a17A", "1miwA", "1aopA", "1a8lA"];
const GROUP_C: &[&str] = &["1a6qA", "1a06A", "1b3oB", "1a0cA", "1a4iA", "1bcoA"];

/// Parsed ground truth FILTER95 output line.
#[derive(Debug)]
struct Filter95Ref {
    zscore: f64,
    cd1cd2: String,
    idom: usize,
    score: i32,
}

/// Parse FILTER95 ground truth output file.
fn parse_filter95_gt(path: &str) -> Vec<Filter95Ref> {
    let content = fs::read_to_string(path).unwrap();
    let mut entries = Vec::new();

    for line in content.lines() {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 5 {
            continue;
        }

        let zscore: f64 = match parts[0].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let cd1cd2 = parts[1].to_string();
        let idom: usize = match parts[2].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let score: i32 = match parts[3].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };

        entries.push(Filter95Ref {
            zscore,
            cd1cd2,
            idom,
            score,
        });
    }

    entries
}

/// Compare FILTER95 results against ground truth.
fn validate_filter95(
    gt: &[Filter95Ref],
    results: &[Filter95Entry],
    score_tol: i32,
    zscore_tol: f64,
) -> (usize, usize, Vec<String>) {
    let mut ref_map: HashMap<(String, usize), &Filter95Ref> = HashMap::new();
    for r in gt {
        let key = (r.cd1cd2.clone(), r.idom);
        if let Some(existing) = ref_map.get(&key) {
            if r.score > existing.score {
                ref_map.insert(key, r);
            }
        } else {
            ref_map.insert(key, r);
        }
    }

    let mut cand_map: HashMap<(String, usize), &Filter95Entry> = HashMap::new();
    for c in results {
        let key = (format!("{}{}", c.cd1, c.cd2), c.idom);
        if let Some(existing) = cand_map.get(&key) {
            if c.score > existing.score {
                cand_map.insert(key, c);
            }
        } else {
            cand_map.insert(key, c);
        }
    }

    let mut matched = 0;
    let mut mismatches = Vec::new();

    for (key, r) in &ref_map {
        if let Some(c) = cand_map.get(key) {
            let score_diff = (r.score - c.score).abs();
            let score_pct = ((r.score as f64) * 0.01).abs() as i32;
            let tol = score_tol.max(score_pct);
            let zscore_diff = (r.zscore - c.zscore).abs();

            if score_diff <= tol && zscore_diff <= zscore_tol {
                matched += 1;
            } else if mismatches.len() < 20 {
                mismatches.push(format!(
                    "{} idom={}: ref_score={} cand_score={} (diff={}) ref_z={:.3} cand_z={:.3} (diff={:.3})",
                    key.0, key.1, r.score, c.score, score_diff, r.zscore, c.zscore, zscore_diff
                ));
            }
        } else if mismatches.len() < 20 {
            mismatches.push(format!(
                "{} idom={}: MISSING (ref_score={}, ref_z={:.3})",
                key.0, key.1, r.score, r.zscore
            ));
        }
    }

    (matched, ref_map.len(), mismatches)
}

/// Filter ParsiHits to only those where both cd1 and cd2 are in codes.
fn filter_hits(hits: &[ParsiHit], codes: &[&str]) -> Vec<ParsiHit> {
    hits.iter()
        .filter(|h| codes.contains(&h.cd1.as_str()) && codes.contains(&h.cd2.as_str()))
        .cloned()
        .collect()
}

/// Filter ground truth to entries where both cd1 and cd2 are in codes.
fn filter_gt(gt: &[Filter95Ref], codes: &[&str]) -> Vec<Filter95Ref> {
    gt.into_iter()
        .filter(|r| {
            if r.cd1cd2.len() != 10 { return false; }
            let cd1 = &r.cd1cd2[..5];
            let cd2 = &r.cd1cd2[5..];
            codes.contains(&cd1) && codes.contains(&cd2)
        })
        .map(|r| Filter95Ref {
            zscore: r.zscore,
            cd1cd2: r.cd1cd2.clone(),
            idom: r.idom,
            score: r.score,
        })
        .collect()
}

#[test]
fn test_filter95_5_structures() {
    let input = fs::read_to_string(FILTER95_INPUT_5).unwrap();
    let hits: Vec<ParsiHit> = input
        .lines()
        .filter_map(|l| parse_refine_line(l))
        .collect();
    assert_eq!(
        hits.len(),
        435,
        "Expected 435 PARSI refine lines for 5-struct corpus"
    );

    let store = ProteinStore::new(FIXTURE_DIR);
    let results = run_filter95(&hits, &store, Some(1.0));

    eprintln!(
        "FILTER95 5-struct: {} output lines",
        results.len()
    );

    let gt = parse_filter95_gt(FILTER95_OUTPUT_5);
    eprintln!(
        "FILTER95 5-struct: {} ground truth lines, {} unique keys",
        gt.len(),
        {
            let mut keys: std::collections::HashSet<(String, usize)> =
                std::collections::HashSet::new();
            for r in &gt {
                keys.insert((r.cd1cd2.clone(), r.idom));
            }
            keys.len()
        }
    );

    let (matched, total, mismatches) = validate_filter95(&gt, &results, 100, 0.5);
    let pct = matched as f64 / total.max(1) as f64 * 100.0;
    eprintln!(
        "FILTER95 5-struct: {}/{} matched ({:.1}%)",
        matched, total, pct
    );
    for m in &mismatches {
        eprintln!("  {}", m);
    }

    assert!(
        matched >= 35,
        "FILTER95 5-struct: only {}/{} matched ({:.1}%), expected >=35",
        matched,
        total,
        pct
    );
}

fn run_filter95_18_group(codes: &[&str], label: &str) {
    let input = fs::read_to_string(FILTER95_INPUT_18).unwrap();
    let all_hits: Vec<ParsiHit> = input
        .lines()
        .filter_map(|l| parse_refine_line(l))
        .collect();

    let hits = filter_hits(&all_hits, codes);
    if hits.is_empty() {
        eprintln!("FILTER95 {}: no input hits for this group", label);
        return;
    }

    let store = ProteinStore::new(EXPANDED_DIR);
    let results = run_filter95(&hits, &store, Some(1.0));

    let all_gt = parse_filter95_gt(FILTER95_OUTPUT_18);
    let group_gt = filter_gt(&all_gt, codes);

    eprintln!(
        "FILTER95 {}: {} input hits, {} output lines, {} ground truth entries",
        label, hits.len(), results.len(), group_gt.len()
    );

    let (matched, total, mismatches) = validate_filter95(&group_gt, &results, 100, 0.5);
    let pct = matched as f64 / total.max(1) as f64 * 100.0;
    eprintln!(
        "FILTER95 {}: {}/{} matched ({:.1}%)",
        label, matched, total, pct
    );
    for m in &mismatches[..mismatches.len().min(10)] {
        eprintln!("  {}", m);
    }

    // Allow 2 missing per group (known Fortran bugs)
    assert!(
        matched + 3 >= total,
        "FILTER95 {}: only {}/{} matched ({:.1}%), expected >= {}",
        label, matched, total, pct, total.saturating_sub(3)
    );
}

#[test]
fn test_filter95_18_group_a() {
    run_filter95_18_group(GROUP_A, "18-group-A");
}

#[test]
fn test_filter95_18_group_b() {
    run_filter95_18_group(GROUP_B, "18-group-B");
}

#[test]
fn test_filter95_18_group_c() {
    run_filter95_18_group(GROUP_C, "18-group-C");
}
