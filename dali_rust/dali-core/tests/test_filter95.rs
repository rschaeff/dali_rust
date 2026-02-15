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

/// Parsed ground truth FILTER95 output line.
#[derive(Debug)]
struct Filter95Ref {
    zscore: f64,
    cd1cd2: String,
    idom: usize,
    score: i32,
    nseg: usize,
}

/// Parse FILTER95 ground truth output file.
///
/// Format: `zscore cd1cd2 idom score nseg ranges...`
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
        let nseg: usize = match parts[4].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };

        entries.push(Filter95Ref {
            zscore,
            cd1cd2,
            idom,
            score,
            nseg,
        });
    }

    entries
}

/// Compare FILTER95 results against ground truth.
///
/// Groups by (cd1cd2, idom), keeps best score per key.
/// Returns (matched, total_ref_keys, mismatches).
fn validate_filter95(
    gt: &[Filter95Ref],
    results: &[Filter95Entry],
    score_tol: i32,
    zscore_tol: f64,
) -> (usize, usize, Vec<String>) {
    // Group gt by (cd1cd2, idom), keep best score
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

    // Group results by (cd1+cd2, idom), keep best score
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

    // Expected: 38/38 unique keys matched
    assert!(
        matched >= 35,
        "FILTER95 5-struct: only {}/{} matched ({:.1}%), expected >=35",
        matched,
        total,
        pct
    );
}

#[test]
fn test_filter95_18_structures() {
    let input = fs::read_to_string(FILTER95_INPUT_18).unwrap();
    let hits: Vec<ParsiHit> = input
        .lines()
        .filter_map(|l| parse_refine_line(l))
        .collect();
    assert_eq!(
        hits.len(),
        13515,
        "Expected 13515 PARSI refine lines for 18-struct corpus"
    );

    let store = ProteinStore::new(EXPANDED_DIR);
    let results = run_filter95(&hits, &store, Some(1.0));

    eprintln!(
        "FILTER95 18-struct: {} output lines",
        results.len()
    );

    let gt = parse_filter95_gt(FILTER95_OUTPUT_18);
    eprintln!(
        "FILTER95 18-struct: {} ground truth lines",
        gt.len()
    );

    let (matched, total, mismatches) = validate_filter95(&gt, &results, 100, 0.5);
    let pct = matched as f64 / total.max(1) as f64 * 100.0;
    eprintln!(
        "FILTER95 18-struct: {}/{} matched ({:.1}%)",
        matched, total, pct
    );
    for m in &mismatches[..mismatches.len().min(10)] {
        eprintln!("  {}", m);
    }

    // Expected: ~293/293 (2 Fortran bug entries correctly rejected)
    assert!(
        matched >= 250,
        "FILTER95 18-struct: only {}/{} matched ({:.1}%), expected >=250",
        matched,
        total,
        pct
    );
}
