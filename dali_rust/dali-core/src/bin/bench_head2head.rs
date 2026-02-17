//! Standalone benchmark binary for head-to-head Rust vs Fortran timing.
//!
//! Runs the Rust pipeline on specified protein pairs and outputs JSON timing.
//! Designed to be called from bench_compare.py alongside the Fortran binary.
//!
//! Usage:
//!   bench_head2head <dat_dir> <method> <cd1> <cd2> [<cd1> <cd2> ...]
//!
//! Methods:
//!   WOLF       - WOLF compare only
//!   PARSI      - PARSI search only
//!   PIPELINE   - Full compare_pair (WOLF+PARSI, both directions)
//!   WOLF_PATH  - Wolf path: WOLF → DP → DALICON → DP
//!   PARSI_PATH - Parsi path: PARSI → FILTER95 → pipe96 → DALICON → DP
//!
//! Output (JSON):
//!   { "method": "...", "pairs": [{ "cd1": "...", "cd2": "...",
//!     "elapsed_ms": ..., "n_results": ... }], "total_ms": ... }

use std::env;
use std::time::Instant;

use dali_core::store::ProteinStore;
use dali_core::{wolf, parsi, pipeline};

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 5 || (args.len() - 3) % 2 != 0 {
        eprintln!("Usage: bench_head2head <dat_dir> <method> <cd1> <cd2> [<cd1> <cd2> ...]");
        eprintln!("Methods: WOLF, PARSI, PIPELINE, WOLF_PATH, PARSI_PATH");
        std::process::exit(1);
    }

    let dat_dir = &args[1];
    let method = &args[2];

    let mut pair_args: Vec<(&str, &str)> = Vec::new();
    let mut i = 3;
    while i + 1 < args.len() {
        pair_args.push((&args[i], &args[i + 1]));
        i += 2;
    }

    let store = ProteinStore::new(dat_dir);

    // Warmup: run first pair once to populate caches and trigger JIT/page faults
    if !pair_args.is_empty() {
        let (cd1, cd2) = pair_args[0];
        match method.as_str() {
            "WOLF" => {
                let cd1_path = format!("{}/{}.dat", dat_dir, cd1);
                let cd2_path = format!("{}/{}.dat", dat_dir, cd2);
                if let Ok(cd1_data) = wolf::setup_protein(&cd1_path) {
                    if let Ok(cd2_data) = wolf::setup_protein(&cd2_path) {
                        let mut grid = wolf::spatial_hash::SpatialHashGrid::new();
                        wolf::load_protein(&mut grid, &cd1_data);
                        let _ = wolf::wolf_compare(&cd1_data, &cd2_data, &grid);
                    }
                }
            }
            "PARSI" => {
                let mut cache = None;
                let _ = parsi::dowork_parsi(cd1, cd2, dat_dir, None, true, &mut cache);
            }
            _ => {
                let _ = pipeline::compare_pair(cd1, cd2, &store);
            }
        }
    }

    let total_start = Instant::now();
    let mut pair_results = Vec::new();

    for (cd1, cd2) in &pair_args {
        let pair_start = Instant::now();
        let n_results;

        match method.as_str() {
            "WOLF" => {
                let cd1_path = format!("{}/{}.dat", dat_dir, cd1);
                let cd2_path = format!("{}/{}.dat", dat_dir, cd2);
                let cd1_data = wolf::setup_protein(&cd1_path).unwrap();
                let cd2_data = wolf::setup_protein(&cd2_path).unwrap();
                let mut grid = wolf::spatial_hash::SpatialHashGrid::new();
                wolf::load_protein(&mut grid, &cd1_data);
                let result = wolf::wolf_compare(&cd1_data, &cd2_data, &grid);
                n_results = if result.is_some() { 1 } else { 0 };
            }
            "PARSI" => {
                let mut cache = None;
                let hits = parsi::dowork_parsi(cd1, cd2, dat_dir, None, true, &mut cache);
                n_results = hits.len();
            }
            "PIPELINE" => {
                let results = pipeline::compare_pair(cd1, cd2, &store);
                n_results = results.len();
            }
            "WOLF_PATH" => {
                let targets: Vec<&str> = vec![cd2];
                let results = pipeline::run_wolf_path(cd1, &targets, &store);
                n_results = results.len();
            }
            "PARSI_PATH" => {
                let targets: Vec<&str> = vec![cd2];
                let results = pipeline::run_parsi_path(cd1, &targets, &store);
                n_results = results.len();
            }
            other => {
                eprintln!("Unknown method: {}", other);
                std::process::exit(1);
            }
        }

        let elapsed = pair_start.elapsed();
        pair_results.push(format!(
            r#"    {{"cd1": "{}", "cd2": "{}", "elapsed_ms": {:.3}, "n_results": {}}}"#,
            cd1, cd2, elapsed.as_secs_f64() * 1000.0, n_results,
        ));
    }

    let total_elapsed = total_start.elapsed();

    // Output JSON manually to avoid serde dependency in hot path
    println!("{{");
    println!(r#"  "method": "{}","#, method);
    println!(r#"  "backend": "rust","#);
    println!("  \"pairs\": [");
    for (i, r) in pair_results.iter().enumerate() {
        if i < pair_results.len() - 1 {
            println!("{},", r);
        } else {
            println!("{}", r);
        }
    }
    println!("  ],");
    println!(r#"  "total_ms": {:.3}"#, total_elapsed.as_secs_f64() * 1000.0);
    println!("}}");
}
