//! Pipeline orchestration module.
//!
//! Wires together:
//!   Path A (Wolf): WOLF → DP(zcut=2.0) → DALICON(lfitz=T) → DP(zcut=2.0)
//!   Path B (Parsi): PARSI → FILTER95(zcut=1.0) → pipe96(top 2) → DALICON(lfitz=T) → DP(zcut=2.0)
//!
//! Both paths run in forward (query=cd1) and reverse (query=cd2) directions.
//! Self-comparisons are filtered from final output.

use std::collections::HashMap;

use crate::types::{AlignmentBlock, DccpEntry};
use crate::store::ProteinStore;
use crate::wolf;
use crate::dp;
use crate::dalicon::{self, DaliconRecord};
use crate::parsi;
use crate::filter95::{self, Filter95Entry};

/// Max entries per (cd1,cd2) key from FILTER95 → DALICON.
const PIPENBEST: usize = 2;

// ---------------------------------------------------------------------------
// Format converters
// ---------------------------------------------------------------------------

/// Convert DccpEntry blocks to DaliconRecord format.
///
/// Flattens AlignmentBlocks into the values array that `convert_prealignment`
/// expects: values[0..2*nblock] = cd1 (l1,r1) pairs, values[2*nblock..4*nblock] = cd2 (l2,r2) pairs.
fn dccp_to_dalicon_record(entry: &DccpEntry) -> DaliconRecord {
    let nblock = entry.blocks.len();
    let mut values = Vec::with_capacity(nblock * 4);
    // cd1 ranges
    for b in &entry.blocks {
        values.push(b.l1 as i32);
        values.push(b.r1 as i32);
    }
    // cd2 ranges
    for b in &entry.blocks {
        values.push(b.l2 as i32);
        values.push(b.r2 as i32);
    }
    DaliconRecord {
        cd1: entry.cd1.clone(),
        cd2: entry.cd2.clone(),
        nblock,
        values,
    }
}

/// Convert Filter95Entry list to DaliconRecords via pipe96 filtering.
///
/// Matches: `sort -nr | uniq | pipe96-free.pl zcut nbest`
/// - Sort by zscore descending
/// - Keep top `nbest` per (cd1+cd2) key
/// - Filter by z-score cutoff
fn filter95_to_dalicon_records(
    entries: &[Filter95Entry],
    zcut: f64,
    nbest: usize,
) -> Vec<DaliconRecord> {
    let mut sorted: Vec<&Filter95Entry> = entries.iter().collect();
    sorted.sort_by(|a, b| b.zscore.partial_cmp(&a.zscore).unwrap_or(std::cmp::Ordering::Equal));

    let mut key_counts: HashMap<String, usize> = HashMap::new();
    let mut records = Vec::new();

    for e in sorted {
        if e.zscore < zcut {
            continue;
        }
        let key = format!("{}{}", e.cd1, e.cd2);
        let count = key_counts.entry(key).or_insert(0);
        *count += 1;
        if *count > nbest {
            continue;
        }

        records.push(DaliconRecord {
            cd1: e.cd1.clone(),
            cd2: e.cd2.clone(),
            nblock: e.nseg,
            values: e.ranges.clone(),
        });
    }

    records
}

// ---------------------------------------------------------------------------
// Pipeline paths
// ---------------------------------------------------------------------------

/// Wolf path: WOLF → DP(zcut=2.0) → DALICON(lfitz=T) → DP(zcut=2.0).
pub fn run_wolf_path(
    query: &str,
    targets: &[&str],
    store: &ProteinStore,
) -> Vec<DccpEntry> {
    let dat_dir = store.dat_dir();

    // Step 1: WOLF — setup query protein and build spatial hash grid
    let cd1_path = format!("{}/{}.dat", dat_dir, query);
    let cd1_data = match wolf::setup_protein(&cd1_path) {
        Ok(d) => d,
        Err(_) => return Vec::new(),
    };

    if cd1_data.protein.nseg <= 2 {
        return Vec::new();
    }

    let mut grid = wolf::spatial_hash::SpatialHashGrid::new();
    wolf::load_protein(&mut grid, &cd1_data);

    let mut wolf_results = Vec::new();
    for &target in targets {
        let cd2_path = format!("{}/{}.dat", dat_dir, target);
        let cd2_data = match wolf::setup_protein(&cd2_path) {
            Ok(d) => d,
            Err(_) => continue,
        };
        if let Some(result) = wolf::wolf_compare(&cd1_data, &cd2_data, &grid) {
            wolf_results.push(result);
        }
    }

    if wolf_results.is_empty() {
        return Vec::new();
    }

    // Step 2: DP on WOLF output
    let alignments: Vec<(String, String, Vec<AlignmentBlock>)> = wolf_results
        .into_iter()
        .map(|r| (r.cd1, r.cd2, r.blocks))
        .collect();

    let dp_results = dp::run_dp(&alignments, store, 2.0);
    if dp_results.is_empty() {
        return Vec::new();
    }

    // Step 3: DALICON on DP output (lfitz=true)
    let records: Vec<DaliconRecord> = dp_results.iter().map(dccp_to_dalicon_record).collect();
    let dalicon_output = dalicon::run_dalicon(&records, store);

    if dalicon_output.is_empty() {
        return Vec::new();
    }

    // Step 4: Final DP
    let final_alignments: Vec<(String, String, Vec<AlignmentBlock>)> = dalicon_output
        .into_iter()
        .map(|(cd1, cd2, _nblock, blocks)| (cd1, cd2, blocks))
        .collect();

    dp::run_dp(&final_alignments, store, 2.0)
}

/// Parsi path: PARSI → FILTER95(zcut=1.0) → pipe96(top 2) → DALICON(lfitz=T) → DP(zcut=2.0).
pub fn run_parsi_path(
    query: &str,
    targets: &[&str],
    store: &ProteinStore,
) -> Vec<DccpEntry> {
    let dat_dir = store.dat_dir();

    // Step 1: PARSI
    let mut all_hits = Vec::new();
    let mut cd1_cache: Option<parsi::ParsiCd1Cache> = None;

    for &target in targets {
        let hits = parsi::dowork_parsi(
            query, target, dat_dir, None, true, &mut cd1_cache,
        );
        all_hits.extend(hits);
    }

    if all_hits.is_empty() {
        return Vec::new();
    }

    // Step 2: FILTER95
    let filter_results = filter95::run_filter95(&all_hits, store, Some(1.0));
    if filter_results.is_empty() {
        return Vec::new();
    }

    // Step 3: pipe96 → DALICON input records
    let records = filter95_to_dalicon_records(&filter_results, 1.0, PIPENBEST);
    if records.is_empty() {
        return Vec::new();
    }

    // Step 4: DALICON (lfitz=true)
    let dalicon_output = dalicon::run_dalicon(&records, store);
    if dalicon_output.is_empty() {
        return Vec::new();
    }

    // Step 5: Final DP
    let final_alignments: Vec<(String, String, Vec<AlignmentBlock>)> = dalicon_output
        .into_iter()
        .map(|(cd1, cd2, _nblock, blocks)| (cd1, cd2, blocks))
        .collect();

    dp::run_dp(&final_alignments, store, 2.0)
}

// ---------------------------------------------------------------------------
// Top-level entry point
// ---------------------------------------------------------------------------

/// Full comparison of two structures.
///
/// Runs wolf and parsi paths in both forward (query=cd1) and reverse
/// (query=cd2) directions. Self-comparisons are filtered out.
pub fn compare_pair(
    cd1: &str,
    cd2: &str,
    store: &ProteinStore,
) -> Vec<DccpEntry> {
    let targets: Vec<&str> = vec![cd1, cd2];
    let mut all_results = Vec::new();

    // Forward direction: query=cd1
    all_results.extend(run_wolf_path(cd1, &targets, store));
    all_results.extend(run_parsi_path(cd1, &targets, store));

    // Reverse direction: query=cd2
    all_results.extend(run_wolf_path(cd2, &targets, store));
    all_results.extend(run_parsi_path(cd2, &targets, store));

    // Filter out self-comparisons
    all_results.retain(|r| r.cd1 != r.cd2);

    all_results
}
