use std::fs;

use dali_core::pipeline::compare_pair;
use dali_core::ProteinStore;
use dali_core::DccpEntry;

const FIXTURE_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../..",
    "/validation/fixtures/structures"
);
const E2E_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../..",
    "/validation/fixtures/e2e"
);

/// Parsed reference DCCP entry from ground truth file.
#[derive(Debug)]
struct DccpRef {
    score: f64,
    zmax: f64,
    cd1: String,
    cd2: String,
}

/// Parse a reference DCCP file into entries (skipping alignment lines).
fn parse_dccp_ref(path: &str) -> Vec<DccpRef> {
    let content = fs::read_to_string(path).unwrap();
    let mut entries = Vec::new();

    for line in content.lines() {
        let trimmed = line.trim();
        if !trimmed.starts_with("DCCP") {
            continue;
        }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() < 10 {
            continue;
        }

        // Format: DCCP 1 score rmsd 0 zmax 0 nblock cd1 cd2
        let score: f64 = match parts[2].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let zmax: f64 = match parts[5].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let cd1 = parts[8].to_string();
        let cd2 = parts[9].to_string();

        entries.push(DccpRef { score, zmax, cd1, cd2 });
    }

    entries
}

/// Compare candidate DccpEntry results against reference, grouped by direction.
///
/// Returns (matched, total, mismatch details).
fn validate_pair(
    candidates: &[DccpEntry],
    refs: &[DccpRef],
    cd1: &str,
    cd2: &str,
    score_tol: f64,
) -> (usize, usize, Vec<String>) {
    let mut matched = 0;
    let mut total = 0;
    let mut details = Vec::new();

    for (dir_cd1, dir_cd2) in [(cd1, cd2), (cd2, cd1)] {
        let ref_dir: Vec<&DccpRef> = refs.iter()
            .filter(|r| r.cd1 == dir_cd1 && r.cd2 == dir_cd2)
            .collect();
        let cand_dir: Vec<&DccpEntry> = candidates.iter()
            .filter(|r| r.cd1 == dir_cd1 && r.cd2 == dir_cd2)
            .collect();

        if ref_dir.is_empty() && cand_dir.is_empty() {
            continue;
        }

        total += ref_dir.len().max(cand_dir.len());

        if ref_dir.is_empty() {
            details.push(format!(
                "  {}→{}: {} EXTRA results (ref has 0)",
                dir_cd1, dir_cd2, cand_dir.len()
            ));
            continue;
        }
        if cand_dir.is_empty() {
            details.push(format!(
                "  {}→{}: MISSING (ref has {})",
                dir_cd1, dir_cd2, ref_dir.len()
            ));
            continue;
        }

        for i in 0..ref_dir.len().min(cand_dir.len()) {
            let r = &ref_dir[i];
            let c = &cand_dir[i];
            let score_diff = (r.score - c.score).abs();

            if score_diff <= score_tol {
                matched += 1;
                details.push(format!(
                    "  {}→{}[{}]: MATCH score={:.1} (ref={:.1} diff={:.1})",
                    dir_cd1, dir_cd2, i, c.score, r.score, score_diff
                ));
            } else {
                details.push(format!(
                    "  {}→{}[{}]: DIFF  score={:.1} (ref={:.1} diff={:.1})",
                    dir_cd1, dir_cd2, i, c.score, r.score, score_diff
                ));
            }
        }
    }

    (matched, total, details)
}

#[test]
fn test_e2e_pipeline() {
    let store = ProteinStore::new(FIXTURE_DIR);
    let score_tol = 2.0;

    let pairs = [
        ("101mA", "1a00A"),
        ("101mA", "1binA"),
        ("1a87A", "1allA"),
    ];

    let mut total_matched = 0;
    let mut total_count = 0;

    for (cd1, cd2) in &pairs {
        let ref_path = format!("{}/{}_vs_{}/{}.dccp", E2E_DIR, cd1, cd2, cd1);
        let refs = parse_dccp_ref(&ref_path);
        // Filter self-comparisons from reference
        let refs: Vec<DccpRef> = refs.into_iter()
            .filter(|r| r.cd1 != r.cd2)
            .collect();

        eprintln!("\n  --- {} vs {} ---", cd1, cd2);
        eprintln!("  Reference: {} DCCP entries", refs.len());

        let results = compare_pair(cd1, cd2, &store);
        eprintln!("  Pipeline:  {} DCCP results", results.len());

        let (matched, count, details) = validate_pair(
            &results, &refs, cd1, cd2, score_tol,
        );

        for d in &details {
            eprintln!("{}", d);
        }

        total_matched += matched;
        total_count += count;
    }

    eprintln!("\nE2E TOTAL: {}/{} matched (score_tol={})", total_matched, total_count, score_tol);

    assert!(
        total_matched >= 8,
        "E2E: only {}/{} matched, expected >=8",
        total_matched, total_count,
    );
}
