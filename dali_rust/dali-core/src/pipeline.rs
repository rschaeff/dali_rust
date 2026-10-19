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
use crate::io::{import_pdb_text, ImportError};
use crate::store::ProteinStore;
use crate::numerics::compute_transform;
use crate::wolf;
use crate::soap;
use crate::dp;
use crate::dalicon::{self, DaliconRecord};
use crate::parsi;
use crate::filter95::{self, Filter95Entry};

/// Max entries per (cd1,cd2) key from FILTER95 → DALICON.
const PIPENBEST_DEFAULT: usize = 2;

/// FILTER95 z-score cutoff, matching `pipe96-free.pl zcut`.
const FILTER95_ZCUT_DEFAULT: f64 = 1.0;

// ---------------------------------------------------------------------------
// Sensitivity tunables
// ---------------------------------------------------------------------------
//
// The three gates below decide how many candidates survive to the next stage. Their
// defaults reproduce the Fortran pipeline exactly and are what every previous benchmark
// ran with; they are overridable only so a run can separate "the alignment was never
// found" from "the alignment was found and then thresholded away", which a fixed binary
// cannot distinguish. Setting DALI_ZCUT_DP=0 and DALI_ZCUT_FILTER95=0 emits every
// candidate that reaches scoring, at the cost of a great deal of noise.
//
//   DALI_ZCUT_DP        DP output cutoff, both passes           (default 2.0)
//   DALI_ZCUT_FILTER95  FILTER95 + pipe96 cutoff                (default 1.0)
//   DALI_PIPENBEST      pipe96 entries kept per (cd1,cd2) key   (default 2)
//
// Read once per process and cached, so changing the variable mid-run does nothing.

/// dali.pl $MINSSE: queries with fewer SSEs go to SOAP.
const MINSSE: usize = 3;

fn env_or<T: std::str::FromStr>(key: &str, default: T) -> T {
    std::env::var(key).ok().and_then(|s| s.parse().ok()).unwrap_or(default)
}

fn zcut_dp() -> f64 {
    static V: std::sync::OnceLock<f64> = std::sync::OnceLock::new();
    *V.get_or_init(|| env_or("DALI_ZCUT_DP", dp::ZCUT_DEFAULT))
}

fn zcut_filter95() -> f64 {
    static V: std::sync::OnceLock<f64> = std::sync::OnceLock::new();
    *V.get_or_init(|| env_or("DALI_ZCUT_FILTER95", FILTER95_ZCUT_DEFAULT))
}

fn pipenbest() -> usize {
    static V: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *V.get_or_init(|| env_or("DALI_PIPENBEST", PIPENBEST_DEFAULT))
}

// ---------------------------------------------------------------------------
// Format converters
// ---------------------------------------------------------------------------

/// Convert DccpEntry blocks to DaliconRecord format.
///
/// Flattens AlignmentBlocks into the values array that `convert_prealignment`
/// expects: values[0..2*nblock] = cd1 (l1,r1) pairs, values[2*nblock..4*nblock] = cd2 (l2,r2) pairs.
pub fn dccp_to_dalicon_record(entry: &DccpEntry) -> DaliconRecord {
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
/// - Sort by zscore descending (ties: the rest of the line, also descending)
/// - `uniq` drops repeated identical lines
/// - Keep top `nbest` per (cd1+cd2) key, with zscore >= zcut
/// - pipe96 strips every `-99` token textually; the remaining tokens are
///   split in half between cd1 and cd2 ranges, nblock = count/4
pub fn filter95_to_dalicon_records(
    entries: &[Filter95Entry],
    zcut: f64,
    nbest: usize,
) -> Vec<DaliconRecord> {
    let mut sorted: Vec<&Filter95Entry> = entries.iter().collect();
    sorted.sort_by(|a, b| {
        b.zscore
            .partial_cmp(&a.zscore)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| (&b.cd1, &b.cd2).cmp(&(&a.cd1, &a.cd2)))
            .then_with(|| b.idom.cmp(&a.idom))
            .then_with(|| b.score.cmp(&a.score))
            .then_with(|| b.nseg.cmp(&a.nseg))
            .then_with(|| b.ranges.cmp(&a.ranges))
    });
    sorted.dedup_by(|b, a| {
        a.zscore == b.zscore
            && a.cd1 == b.cd1
            && a.cd2 == b.cd2
            && a.idom == b.idom
            && a.score == b.score
            && a.nseg == b.nseg
            && a.ranges == b.ranges
    });

    let mut key_counts: HashMap<String, usize> = HashMap::new();
    let mut records = Vec::new();

    for e in sorted {
        let key = format!("{}{}", e.cd1, e.cd2);
        let count = key_counts.entry(key).or_insert(0);
        *count += 1;
        if *count > nbest || e.zscore < zcut {
            continue;
        }

        let kept: Vec<i32> = e.ranges.iter().copied().filter(|&v| v != -99).collect();
        let n = kept.len() / 4;
        let mut values = kept[..2 * n].to_vec();
        values.extend_from_slice(&kept[2 * n..4 * n]);
        records.push(DaliconRecord {
            cd1: e.cd1.clone(),
            cd2: e.cd2.clone(),
            nblock: n,
            values,
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

    wolfitz_to_dccp(wolf_results, store)
}

/// SOAP path (query nsse < 3): SOAP → DP(zcut=2.0) → DALICON(lfitz=T) → DP(zcut=2.0).
///
/// As serialcompare, the query's CA coordinates are loaded once and stay
/// transformed from one target to the next.
pub fn run_soap_path(
    query: &str,
    targets: &[&str],
    store: &ProteinStore,
) -> Vec<DccpEntry> {
    let qprot = match store.get_protein(query) {
        Ok(p) => p,
        Err(_) => return Vec::new(),
    };
    let mut q = soap::SoapQuery::new(&qprot);
    let mut results = Vec::new();
    for &target in targets {
        let tprot = match store.get_protein(target) {
            Ok(p) => p,
            Err(_) => continue,
        };
        if let Some(r) = soap::soap_compare(&mut q, &tprot) {
            results.push(r);
        }
    }
    wolfitz_to_dccp(results, store)
}

/// What mpidali.pm's `wolf` does with WOLFITZ lines from WOLF or SOAP:
/// DP → DALICON(lfitz=T) → DP.
fn wolfitz_to_dccp(wolf_results: Vec<wolf::WolfResult>, store: &ProteinStore) -> Vec<DccpEntry> {
    if wolf_results.is_empty() {
        return Vec::new();
    }

    // Step 2: DP on WOLF output
    let alignments: Vec<(String, String, Vec<AlignmentBlock>)> = wolf_results
        .into_iter()
        .map(|r| (r.cd1, r.cd2, r.blocks))
        .collect();

    let dp_results = dp::run_dp(&alignments, store, zcut_dp());
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

    dp::run_dp(&final_alignments, store, zcut_dp())
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
    let filter_results = filter95::run_filter95(&all_hits, store, Some(zcut_filter95()));
    if filter_results.is_empty() {
        return Vec::new();
    }

    // Step 3: pipe96 → DALICON input records
    let records = filter95_to_dalicon_records(&filter_results, zcut_filter95(), pipenbest());
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

    dp::run_dp(&final_alignments, store, zcut_dp())
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

/// Re-express an alignment from the other side: cd1 and cd2 swap, every block
/// swaps its ranges, and blocks are re-sorted by their new cd1 position.
pub fn swap_orientation(entry: &DccpEntry) -> DccpEntry {
    let mut blocks: Vec<AlignmentBlock> = entry.blocks.iter()
        .map(|b| AlignmentBlock { l1: b.l2, r1: b.r2, l2: b.l1, r2: b.r1 })
        .collect();
    blocks.sort_by_key(|b| b.l1);
    DccpEntry {
        cd1: entry.cd2.clone(),
        cd2: entry.cd1.clone(),
        score: entry.score,
        zscore: entry.zscore,
        rmsd: entry.rmsd,
        blocks,
    }
}

/// One-to-many structural search with query-side amortization.
///
/// Runs WOLF and/or PARSI paths query→targets. With `bothways`, also runs each
/// target as cd1 against the query and swaps those hits into query orientation,
/// which is what `dali.pl` pairwise mode does (do_one_pairwise without --oneway):
/// WOLF and PARSI are both asymmetric, and for a small template inside a large
/// query the template→query direction is often the one that finds it.
///
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
    bothways: bool,
) -> Vec<DccpEntry> {
    let mut all_results = Vec::new();

    // dali.pl do_one_pairwise: a query with nsse < MINSSE (3) runs SOAP
    // only, in the forward direction.
    let soap_query = store.get_protein(query).map_or(false, |p| p.nseg < MINSSE);
    if soap_query {
        all_results.extend(run_soap_path(query, targets, store));
    } else if !skip_wolf {
        all_results.extend(run_wolf_path(query, targets, store));
    }
    if !soap_query {
        all_results.extend(run_parsi_path(query, targets, store));
    }

    if bothways && !soap_query {
        let q = [query];
        for &target in targets.iter().filter(|&&t| t != query) {
            let mut rev = Vec::new();
            if !skip_wolf {
                rev.extend(run_wolf_path(target, &q, store));
            }
            rev.extend(run_parsi_path(target, &q, store));
            all_results.extend(
                rev.iter().filter(|r| r.cd1 == target && r.cd2 == query).map(swap_orientation),
            );
        }
    }

    // Filter self-comparisons
    all_results.retain(|r| r.cd1 != r.cd2);

    // Filter by z-score
    all_results.retain(|r| r.zscore >= z_cut);

    // Deduplicate as dali.pl's report does: sortdccp.pl sorts the DCCP blocks
    // of cd1.dccp on the Z-score *as printed* (f8.1) with Perl's stable sort,
    // and fsspfilter.pl keeps the first block per cd2. all_results is in
    // cd1.dccp order (WOLF, PARSI, then the reverse direction's WOLF, PARSI),
    // so among blocks whose printed Z ties, the earliest wins.
    let printed_z = |z: f64| -> f64 { format!("{:.1}", z).parse().unwrap_or(z) };
    let mut order: Vec<(String, String)> = Vec::new();
    let mut best_map: HashMap<(String, String), DccpEntry> = HashMap::new();
    for entry in all_results {
        let key = (entry.cd1.clone(), entry.cd2.clone());
        let replace = match best_map.get(&key) {
            Some(existing) => printed_z(entry.zscore) > printed_z(existing.zscore),
            None => {
                order.push(key.clone());
                true
            }
        };
        if replace {
            best_map.insert(key, entry);
        }
    }

    let mut results: Vec<DccpEntry> = order.into_iter()
        .filter_map(|k| best_map.remove(&k))
        .collect();
    results.sort_by(|a, b| printed_z(b.zscore).partial_cmp(&printed_z(a.zscore))
        .unwrap_or(std::cmp::Ordering::Equal));
    results
}

// ---------------------------------------------------------------------------
// Iterative search (search + mask loop)
// ---------------------------------------------------------------------------

/// Iterative multi-domain detection: search → take best → mask → repeat.
///
/// Each round runs `search_database` (in both directions if `bothways`), takes the top hit, computes the
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
    bothways: bool,
) -> Vec<SearchHit> {
    let mut hits = Vec::new();
    let mut current_query = query.to_string();

    for round in 0..max_rounds {
        let target_refs: Vec<&str> = targets.iter().copied().collect();
        let results = search_database(&current_query, &target_refs, store, min_zscore, skip_wolf, bothways);

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
            query_resids: alignments.iter().map(|&(q, _)| query_prot.resid_map[q - 1]).collect(),
            qlen: query_prot.nres,
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

// ---------------------------------------------------------------------------
// Iterative search on PDB text (DPAM step 7, Fortran semantics)
// ---------------------------------------------------------------------------

/// Parameters for [`iterative_search_pdb`]; defaults are DPAM step 7's.
#[derive(Debug, Clone)]
pub struct PdbIterParams {
    pub min_aligned: usize,
    pub min_zscore: f64,
    pub min_remaining: usize,
    pub max_rounds: usize,
    pub bothways: bool,
}

impl Default for PdbIterParams {
    fn default() -> Self {
        PdbIterParams { min_aligned: 20, min_zscore: 2.0, min_remaining: 20, max_rounds: usize::MAX, bothways: true }
    }
}

/// DPAM `get_domain_range`: sort, split where the next residue number exceeds the
/// previous by more than max(5, 0.05 * n); returns inclusive (first, last) ranges.
pub fn dpam_domain_ranges(resids: &[i32]) -> Vec<(i32, i32)> {
    let mut r = resids.to_vec();
    r.sort_unstable();
    let cutoff = 5.0f64.max(r.len() as f64 * 0.05);
    let mut segs: Vec<(i32, i32)> = Vec::new();
    for &x in &r {
        match segs.last_mut() {
            Some(s) if (x as f64) <= s.1 as f64 + cutoff => s.1 = x,
            _ => segs.push((x, x)),
        }
    }
    segs
}

/// Residue numbers (columns 23-26 as an integer) of the ATOM records, as DPAM's
/// step 7 reads them; insertion codes and chains are ignored, as there.
fn atom_resids(text: &str) -> Vec<i32> {
    let mut v: Vec<i32> = text.lines()
        .filter(|l| l.starts_with("ATOM"))
        .filter_map(|l| l.get(22..26).and_then(|s| s.trim().parse().ok()))
        .collect();
    v.sort_unstable();
    v.dedup();
    v
}

/// Iterative domain detection exactly as DPAM step 7 drives Fortran DaliLite:
/// every round re-imports the CURRENT query PDB (so DSSP, PUU and the .dat are
/// recomputed for what is left), aligns it against the targets, takes the best hit,
/// removes the residues spanned by `dpam_domain_ranges` of its aligned query residue
/// numbers, and writes the next round's query from the ATOM records that remain.
///
/// Stops when no hit reaches `min_zscore`, the best hit aligns fewer than
/// `min_aligned` residues, or fewer than `min_remaining` residues would remain.
/// Round r's query is stored as `{code}_r{r}`.
pub fn iterative_search_pdb(
    query_pdb: &str,
    chain: &str,
    code: &str,
    targets: &[&str],
    store: &ProteinStore,
    p: &PdbIterParams,
) -> Result<Vec<SearchHit>, ImportError> {
    let mut hits = Vec::new();
    let mut text = query_pdb.to_string();
    for round in 0..p.max_rounds {
        let mut prot = match import_pdb_text(&text, chain, code) {
            Ok(x) => x,
            Err(e) if round == 0 => return Err(e),
            Err(_) => break, // what remains no longer imports (e.g. no complete residue)
        };
        let qcode = format!("{}_r{}", code, round);
        prot.code = qcode.clone();
        store.add_protein(prot)
            .map_err(|e| ImportError::Pdb(crate::io::PdbReadError::Io(e.to_string())))?;
        let qprot = match store.get_protein(&qcode) {
            Ok(x) => x,
            Err(_) => break,
        };

        let results = search_database(&qcode, targets, store, p.min_zscore, false, p.bothways);
        let best = match results.first() {
            Some(b) => b,
            None => break,
        };
        let n_aligned: usize = best.blocks.iter().map(|b| (b.r1 - b.l1 + 1) as usize).sum();
        if n_aligned < p.min_aligned {
            break;
        }
        let template = match store.get_protein(&best.cd2) {
            Ok(t) => t,
            Err(_) => break,
        };
        let (rotation, translation) = match compute_transform(&qprot.ca, &template.ca, &best.blocks) {
            Some(x) => x,
            None => break,
        };
        let mut alignments = Vec::new();
        for b in &best.blocks {
            for k in 0..(b.r1 - b.l1 + 1) as usize {
                let (qi, ti) = (b.l1 as usize - 1 + k, b.l2 as usize - 1 + k);
                if qi < qprot.nres && ti < template.nres {
                    alignments.push((qi + 1, ti + 1));
                }
            }
        }
        let query_resids: Vec<i32> = alignments.iter().map(|&(q, _)| qprot.resid_map[q - 1]).collect();
        let present = atom_resids(&text);
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
            query_resids: query_resids.clone(),
            qlen: present.len(),
        });

        let ranges = dpam_domain_ranges(&query_resids);
        let removed = |r: i32| ranges.iter().any(|&(a, b)| r >= a && r <= b);
        let remain: Vec<i32> = present.iter().copied().filter(|&r| !removed(r)).collect();
        if remain.len() < p.min_remaining {
            break;
        }
        let keep: HashSet<i32> = remain.into_iter().collect();
        text = text.lines()
            .filter(|l| l.starts_with("ATOM"))
            .filter(|l| l.get(22..26).and_then(|s| s.trim().parse::<i32>().ok()).map_or(false, |r| keep.contains(&r)))
            .map(|l| format!("{}\n", l))
            .collect();
    }
    Ok(hits)
}
