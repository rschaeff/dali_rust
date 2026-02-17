//! Criterion benchmarks for dali-core kernels.
//!
//! Run: cargo bench --package dali-core
//!
//! Benchmarks are organized by kernel, each tested across protein sizes
//! from the validation corpus.

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

use dali_core::io::dat;
use dali_core::numerics::{kabsch, scoring, fitz};
use dali_core::store::ProteinStore;
use dali_core::{wolf, dp, dalicon, parsi, filter95, pipeline};
use dali_core::dalicon::DaliconRecord;

// ---------------------------------------------------------------------------
// Fixture paths
// ---------------------------------------------------------------------------

const DAT_DIR_5: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../validation/fixtures/structures"
);
const DAT_DIR_18: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../validation/corpus_expansion/ground_truth_18/structures"
);

/// Representative protein pairs: (cd1, cd2, dat_dir, label).
/// Chosen to span size range: small (62+93=155 res), medium (141+154=295 res),
/// large (297+456=753 res).
fn pairs() -> Vec<(&'static str, &'static str, &'static str, &'static str)> {
    vec![
        ("101mA", "1a00A", DAT_DIR_5, "small_globin"),
        ("1a87A", "1allA", DAT_DIR_5, "medium_globin"),
        ("1aopA", "1a0cA", DAT_DIR_18, "large_ecod"),
    ]
}

// ---------------------------------------------------------------------------
// Numerics: u3b (Kabsch SVD)
// ---------------------------------------------------------------------------

fn bench_u3b(c: &mut Criterion) {
    let mut group = c.benchmark_group("u3b");

    for (cd1, cd2, dat_dir, label) in pairs() {
        let p1 = dat::read_dat(format!("{}/{}.dat", dat_dir, cd1)).unwrap();
        let p2 = dat::read_dat(format!("{}/{}.dat", dat_dir, cd2)).unwrap();
        let n = p1.nres.min(p2.nres);
        let w: Vec<f64> = vec![1.0; n];

        group.bench_with_input(
            BenchmarkId::new("u3b", label),
            &n,
            |b, &n| {
                b.iter(|| {
                    kabsch::u3b(&w, &p1.ca, &p2.ca, n, 1)
                })
            },
        );
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// Numerics: dpgetdist (distance matrix)
// ---------------------------------------------------------------------------

fn bench_dpgetdist(c: &mut Criterion) {
    let mut group = c.benchmark_group("dpgetdist");

    let proteins = [
        ("1aiwA", DAT_DIR_18, "62res"),
        ("101mA", DAT_DIR_5, "154res"),
        ("1a87A", DAT_DIR_5, "297res"),
        ("1aopA", DAT_DIR_18, "456res"),
    ];

    for (code, dat_dir, label) in proteins {
        let p = dat::read_dat(format!("{}/{}.dat", dat_dir, code)).unwrap();
        group.bench_with_input(
            BenchmarkId::new("dpgetdist", label),
            &p.nres,
            |b, _| {
                b.iter(|| scoring::dpgetdist(&p.ca, p.nres))
            },
        );
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// Numerics: fitz (iterative superposition)
// ---------------------------------------------------------------------------

fn bench_fitz(c: &mut Criterion) {
    let mut group = c.benchmark_group("fitz");

    for (cd1, cd2, dat_dir, label) in pairs() {
        let p1 = dat::read_dat(format!("{}/{}.dat", dat_dir, cd1)).unwrap();
        let p2 = dat::read_dat(format!("{}/{}.dat", dat_dir, cd2)).unwrap();

        group.bench_with_input(
            BenchmarkId::new("fitz", label),
            &label,
            |b, _| {
                let mut x = p1.ca.clone();
                b.iter(|| {
                    x.assign(&p1.ca);
                    fitz::fitz(&mut x, &p2.ca, 4.0, 20)
                })
            },
        );
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// WOLF: setup_protein + wolf_compare
// ---------------------------------------------------------------------------

fn bench_wolf(c: &mut Criterion) {
    let mut group = c.benchmark_group("wolf");
    group.sample_size(20);

    for (cd1, cd2, dat_dir, label) in pairs() {
        let cd1_path = format!("{}/{}.dat", dat_dir, cd1);
        let cd2_path = format!("{}/{}.dat", dat_dir, cd2);

        let cd1_data = wolf::setup_protein(&cd1_path).unwrap();
        let cd2_data = wolf::setup_protein(&cd2_path).unwrap();
        let mut grid = wolf::spatial_hash::SpatialHashGrid::new();
        wolf::load_protein(&mut grid, &cd1_data);

        group.bench_with_input(
            BenchmarkId::new("wolf_compare", label),
            &label,
            |b, _| {
                b.iter(|| wolf::wolf_compare(&cd1_data, &cd2_data, &grid))
            },
        );
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// DP: run_dp on WOLF output
// ---------------------------------------------------------------------------

fn bench_dp(c: &mut Criterion) {
    let mut group = c.benchmark_group("dp");
    group.sample_size(20);

    for (cd1, cd2, dat_dir, label) in pairs() {
        let store = ProteinStore::new(dat_dir);

        // Get WOLF output to feed to DP
        let cd1_path = format!("{}/{}.dat", dat_dir, cd1);
        let cd2_path = format!("{}/{}.dat", dat_dir, cd2);
        let cd1_data = wolf::setup_protein(&cd1_path).unwrap();
        let cd2_data = wolf::setup_protein(&cd2_path).unwrap();
        let mut grid = wolf::spatial_hash::SpatialHashGrid::new();
        wolf::load_protein(&mut grid, &cd1_data);

        if let Some(result) = wolf::wolf_compare(&cd1_data, &cd2_data, &grid) {
            let alignments = vec![(result.cd1, result.cd2, result.blocks)];

            group.bench_with_input(
                BenchmarkId::new("run_dp", label),
                &label,
                |b, _| {
                    b.iter(|| dp::run_dp(&alignments, &store, 2.0))
                },
            );
        }
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// DALICON: run_dalicon on DP output
// ---------------------------------------------------------------------------

fn bench_dalicon(c: &mut Criterion) {
    let mut group = c.benchmark_group("dalicon");
    group.sample_size(10);

    for (cd1, cd2, dat_dir, label) in pairs() {
        let store = ProteinStore::new(dat_dir);

        // Run WOLF+DP to get DALICON input
        let cd1_path = format!("{}/{}.dat", dat_dir, cd1);
        let cd2_path = format!("{}/{}.dat", dat_dir, cd2);
        let cd1_data = wolf::setup_protein(&cd1_path).unwrap();
        let cd2_data = wolf::setup_protein(&cd2_path).unwrap();
        let mut grid = wolf::spatial_hash::SpatialHashGrid::new();
        wolf::load_protein(&mut grid, &cd1_data);

        if let Some(wolf_result) = wolf::wolf_compare(&cd1_data, &cd2_data, &grid) {
            let alignments = vec![(wolf_result.cd1, wolf_result.cd2, wolf_result.blocks)];
            let dp_results = dp::run_dp(&alignments, &store, 2.0);

            if !dp_results.is_empty() {
                let records: Vec<DaliconRecord> = dp_results.iter().map(|e| {
                    let nblock = e.blocks.len();
                    let mut values = Vec::with_capacity(nblock * 4);
                    for blk in &e.blocks {
                        values.push(blk.l1 as i32);
                        values.push(blk.r1 as i32);
                    }
                    for blk in &e.blocks {
                        values.push(blk.l2 as i32);
                        values.push(blk.r2 as i32);
                    }
                    DaliconRecord {
                        cd1: e.cd1.clone(),
                        cd2: e.cd2.clone(),
                        nblock,
                        values,
                    }
                }).collect();

                group.bench_with_input(
                    BenchmarkId::new("run_dalicon", label),
                    &label,
                    |b, _| {
                        b.iter(|| dalicon::run_dalicon(&records, &store))
                    },
                );
            }
        }
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// PARSI: dowork_parsi single pair
// ---------------------------------------------------------------------------

fn bench_parsi(c: &mut Criterion) {
    let mut group = c.benchmark_group("parsi");
    group.sample_size(10);

    for (cd1, cd2, dat_dir, label) in pairs() {
        group.bench_with_input(
            BenchmarkId::new("dowork_parsi", label),
            &label,
            |b, _| {
                b.iter(|| {
                    let mut cache = None;
                    parsi::dowork_parsi(cd1, cd2, dat_dir, None, true, &mut cache)
                })
            },
        );
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// FILTER95: run_filter95 on PARSI output
// ---------------------------------------------------------------------------

fn bench_filter95(c: &mut Criterion) {
    let mut group = c.benchmark_group("filter95");
    group.sample_size(20);

    for (cd1, cd2, dat_dir, label) in pairs() {
        let store = ProteinStore::new(dat_dir);
        let mut cache = None;
        let hits = parsi::dowork_parsi(cd1, cd2, dat_dir, None, true, &mut cache);

        if !hits.is_empty() {
            group.bench_with_input(
                BenchmarkId::new("run_filter95", label),
                &label,
                |b, _| {
                    b.iter(|| filter95::run_filter95(&hits, &store, Some(1.0)))
                },
            );
        }
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// Full pipeline: compare_pair
// ---------------------------------------------------------------------------

fn bench_pipeline(c: &mut Criterion) {
    let mut group = c.benchmark_group("pipeline");
    group.sample_size(10);

    for (cd1, cd2, dat_dir, label) in pairs() {
        let store = ProteinStore::new(dat_dir);

        group.bench_with_input(
            BenchmarkId::new("compare_pair", label),
            &label,
            |b, _| {
                b.iter(|| pipeline::compare_pair(cd1, cd2, &store))
            },
        );

        group.bench_with_input(
            BenchmarkId::new("compare_pair_serial", label),
            &label,
            |b, _| {
                b.iter(|| pipeline::compare_pair_serial(cd1, cd2, &store))
            },
        );
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// Pipeline sub-paths: wolf_path and parsi_path
// ---------------------------------------------------------------------------

fn bench_wolf_path(c: &mut Criterion) {
    let mut group = c.benchmark_group("wolf_path");
    group.sample_size(10);

    for (cd1, cd2, dat_dir, label) in pairs() {
        let store = ProteinStore::new(dat_dir);
        let targets: Vec<&str> = vec![cd2];

        group.bench_with_input(
            BenchmarkId::new("run_wolf_path", label),
            &label,
            |b, _| {
                b.iter(|| pipeline::run_wolf_path(cd1, &targets, &store))
            },
        );
    }
    group.finish();
}

fn bench_parsi_path(c: &mut Criterion) {
    let mut group = c.benchmark_group("parsi_path");
    group.sample_size(10);

    for (cd1, cd2, dat_dir, label) in pairs() {
        let store = ProteinStore::new(dat_dir);
        let targets: Vec<&str> = vec![cd2];

        group.bench_with_input(
            BenchmarkId::new("run_parsi_path", label),
            &label,
            |b, _| {
                b.iter(|| pipeline::run_parsi_path(cd1, &targets, &store))
            },
        );
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_u3b,
    bench_dpgetdist,
    bench_fitz,
    bench_wolf,
    bench_dp,
    bench_dalicon,
    bench_parsi,
    bench_filter95,
    bench_pipeline,
    bench_wolf_path,
    bench_parsi_path,
);
criterion_main!(benches);
