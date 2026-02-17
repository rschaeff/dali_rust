//! Pipeline orchestration module.
//!
//! Wires together:
//!   Path A (Wolf): WOLF → DP(zcut=2.0) → DALICON(lfitz=T) → DP(zcut=2.0)
//!   Path B (Parsi): PARSI → FILTER95(zcut=1.0) → pipe96(top 2) → DALICON(lfitz=T) → DP(zcut=2.0)
//!
//! Both paths run in forward (query=cd1) and reverse (query=cd2) directions.
//! Self-comparisons are filtered from final output.

use std::collections::{HashMap, HashSet};

use rayon::prelude::*;

use crate::types::{AlignmentBlock, DccpEntry, SearchHit, mask_protein};
use crate::store::ProteinStore;
use crate::numerics::compute_transform;
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

/// Full comparison of two structures (serial).
///
/// Runs wolf and parsi paths in both forward (query=cd1) and reverse
/// (query=cd2) directions. Self-comparisons are filtered out.
pub fn compare_pair_serial(
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

/// Full comparison of two structures (parallel via Rayon).
///
/// Runs 4 independent paths concurrently:
///   1. Forward wolf  (cd1 → [cd1,cd2])
///   2. Forward parsi (cd1 → [cd1,cd2])
///   3. Reverse wolf  (cd2 → [cd1,cd2])
///   4. Reverse parsi (cd2 → [cd1,cd2])
///
/// Each path has its own mutable state (DaliconCd1State, ParsiCd1Cache).
/// Only the ProteinStore is shared (thread-safe via RwLock + Arc).
pub fn compare_pair(
    cd1: &str,
    cd2: &str,
    store: &ProteinStore,
) -> Vec<DccpEntry> {
    let targets: Vec<&str> = vec![cd1, cd2];

    // Define the 4 independent tasks as (query, path_type) pairs
    let tasks: Vec<(&str, bool)> = vec![
        (cd1, true),   // forward wolf
        (cd1, false),  // forward parsi
        (cd2, true),   // reverse wolf
        (cd2, false),  // reverse parsi
    ];

    let mut all_results: Vec<DccpEntry> = tasks
        .into_par_iter()
        .flat_map(|(query, is_wolf)| {
            if is_wolf {
                run_wolf_path(query, &targets, store)
            } else {
                run_parsi_path(query, &targets, store)
            }
        })
        .collect();

    // Filter out self-comparisons
    all_results.retain(|r| r.cd1 != r.cd2);

    all_results
}

// ---------------------------------------------------------------------------
// Database search (one-to-many)
// ---------------------------------------------------------------------------

/// One-to-many structural search with query-side amortization.
///
/// Runs WOLF and/or PARSI paths (forward only: query→targets).
/// Filters self-comparisons, applies z-score cutoff, deduplicates by
/// keeping the highest z-score per (cd1, cd2) pair, and sorts descending.
///
/// When `skip_wolf` is true, only the PARSI path runs — useful when an
/// upstream prefilter (e.g. FoldSeek) has already narrowed candidates.
pub fn search_database(
    query: &str,
    targets: &[&str],
    store: &ProteinStore,
    z_cut: f64,
    skip_wolf: bool,
) -> Vec<DccpEntry> {
    let mut all_results = Vec::new();

    if !skip_wolf {
        all_results.extend(run_wolf_path(query, targets, store));
    }
    all_results.extend(run_parsi_path(query, targets, store));

    // Filter self-comparisons
    all_results.retain(|r| r.cd1 != r.cd2);

    // Filter by z-score
    all_results.retain(|r| r.zscore >= z_cut);

    // Deduplicate: keep highest z-score per (cd1, cd2) pair
    let mut best_map: HashMap<(String, String), DccpEntry> = HashMap::new();
    for entry in all_results {
        let key = (entry.cd1.clone(), entry.cd2.clone());
        let replace = match best_map.get(&key) {
            Some(existing) => entry.zscore > existing.zscore,
            None => true,
        };
        if replace {
            best_map.insert(key, entry);
        }
    }

    let mut results: Vec<DccpEntry> = best_map.into_values().collect();
    results.sort_by(|a, b| b.zscore.partial_cmp(&a.zscore).unwrap_or(std::cmp::Ordering::Equal));
    results
}

// ---------------------------------------------------------------------------
// Iterative search (search + mask loop)
// ---------------------------------------------------------------------------

/// Iterative multi-domain detection: search → take best → mask → repeat.
///
/// Each round runs `search_database`, takes the top hit, computes the
/// structural transform, masks aligned residues from the query, and repeats
/// until no hit exceeds thresholds or the query becomes too small.
pub fn iterative_search(
    query: &str,
    targets: &[&str],
    store: &ProteinStore,
    min_aligned: usize,
    min_zscore: f64,
    gap_tolerance: usize,
    max_rounds: usize,
    skip_wolf: bool,
) -> Vec<SearchHit> {
    let mut hits = Vec::new();
    let mut current_query = query.to_string();

    for round in 0..max_rounds {
        let target_refs: Vec<&str> = targets.iter().copied().collect();
        let results = search_database(&current_query, &target_refs, store, min_zscore, skip_wolf);

        if results.is_empty() {
            break;
        }

        let best = &results[0];

        // Count aligned residues
        let n_aligned: usize = best.blocks.iter()
            .map(|b| (b.r1 - b.l1 + 1) as usize)
            .sum();
        if n_aligned < min_aligned {
            break;
        }

        // Load proteins for transform computation
        let query_prot = match store.get_protein(&current_query) {
            Ok(p) => p,
            Err(_) => break,
        };
        let template_prot = match store.get_protein(&best.cd2) {
            Ok(p) => p,
            Err(_) => break,
        };

        // Compute rotation + translation
        let (rotation, translation) = match compute_transform(&query_prot.ca, &template_prot.ca, &best.blocks) {
            Some((u, t)) => (u, t),
            None => break,
        };

        // Expand blocks to residue-level alignment pairs (1-based)
        let mut alignments: Vec<(usize, usize)> = Vec::new();
        let mut aligned_indices: HashSet<usize> = HashSet::new(); // 0-based query indices

        for b in &best.blocks {
            let len = (b.r1 - b.l1 + 1) as usize;
            for k in 0..len {
                let q_idx = (b.l1 as usize - 1) + k; // 0-based
                let t_idx = (b.l2 as usize - 1) + k;
                if q_idx < query_prot.nres && t_idx < template_prot.nres {
                    alignments.push((q_idx + 1, t_idx + 1)); // 1-based
                    aligned_indices.insert(q_idx);
                }
            }
        }

        hits.push(SearchHit {
            cd2: best.cd2.clone(),
            zscore: best.zscore,
            score: best.score,
            rmsd: best.rmsd,
            nblock: best.blocks.len(),
            blocks: best.blocks.clone(),
            rotation,
            translation,
            alignments,
            round,
        });

        // Compute mask set: aligned indices + gap-bridged neighbors
        let effective_gap = gap_tolerance.max(n_aligned * 5 / 100);
        let mut mask_set = aligned_indices.clone();

        // Bridge small gaps between aligned regions
        if !aligned_indices.is_empty() {
            let mut sorted_aligned: Vec<usize> = aligned_indices.iter().copied().collect();
            sorted_aligned.sort();
            for window in sorted_aligned.windows(2) {
                let gap = window[1] - window[0];
                if gap > 1 && gap <= effective_gap + 1 {
                    for idx in (window[0] + 1)..window[1] {
                        mask_set.insert(idx);
                    }
                }
            }
        }

        // Compute keep_indices: all residues NOT in mask
        let keep_indices: Vec<usize> = (0..query_prot.nres)
            .filter(|i| !mask_set.contains(i))
            .collect();

        if keep_indices.len() < min_aligned {
            break;
        }

        // Mask and store
        let new_code = format!("{}_r{}", query, round);
        let masked = mask_protein(&query_prot, &keep_indices, &new_code);

        if masked.nseg < 3 {
            // Too few SSEs for WOLF/PARSI to work meaningfully
            break;
        }

        if store.add_protein(masked).is_err() {
            break;
        }
        current_query = new_code;
    }

    hits
}
