# Pre-DPAM Integration: Current State of Rust DALI

**Date:** 2026-02-14
**Purpose:** Document the state of the Rust DALI reimplementation before integrating
it into DPAM v2's Step 7 (Iterative DALI) pipeline.

## 1. What Exists

### Rust Library (`dali-core`)

A fully standalone Rust implementation of the DaliLite.v5 structural comparison
pipeline. It reads PDB/CIF files directly, computes DSSP secondary structure,
builds domain decompositions, and runs the full DALI comparison pipeline with no
external dependencies beyond standard crates.

**Source:** 9,015 lines across 32 `.rs` files in `dali_rust/dali-core/src/`.
**Tests:** 2,075 lines across 10 test files, 60 tests passing.

### Python Bindings (`dali-python`)

PyO3 bindings exposing the Rust library as a Python `dali` module.

**Source:** 213 lines in `dali_rust/dali-python/src/lib.rs`.
**Built with:** maturin 1.12, pyo3 0.23, numpy 0.23.

### Python Types Exposed

```python
dali.Protein         # .code, .nres, .nseg, .na, .nb, .sequence, .ca_coords()
dali.ProteinStore    # __init__(dat_dir), .get_protein(code), __contains__, __len__
dali.DccpEntry       # .cd1, .cd2, .score, .zscore, .rmsd, .nblock, .blocks
dali.AlignmentBlock  # .l1, .r1, .l2, .r2  (1-based residue ranges)
```

### Python Functions Exposed

```python
dali.compare_pair(cd1, cd2, store) -> list[DccpEntry]    # Full pipeline
dali.run_wolf_path(query, targets, store) -> list[DccpEntry]  # WOLF+DP+DALICON
dali.run_parsi_path(query, targets, store) -> list[DccpEntry] # PARSI+FILTER95+DALICON+DP
dali.read_dat(path) -> Protein                           # Load .dat file
dali.import_pdb(path, chain, pdb_code) -> Protein        # PDB/CIF -> Protein
```

## 2. Validation Status

### Module-Level Validation (Rust vs Fortran Ground Truth)

| Module   | 5-struct | 18-struct | Total | Notes |
|----------|----------|-----------|-------|-------|
| WOLF     | 21/21    | 278/279   | 299/300 (99.7%) | 1 float edge case at fitz cutoff |
| DP       | 10/10    | 24/24     | 34/34 (100%) | |
| DALICON  | 10/10    | 24/24     | 34/34 (100%) | |
| PARSI    | 379/435  | 9616/13515| 9995/13950 (71.6%) | Search divergence, not bugs |
| FILTER95 | 38/38    | 292/295   | 330/333 (99.1%) | 2 Fortran bugs correctly rejected |
| Import   | 42/42    | 279/303   | 321/345 (93.0%) | DSSP reimpl variations |
| **E2E**  | **10/10**|           | **10/10 (100%)** | 0.0 score difference |

### End-to-End Validation

Three reference pairs tested through the full pipeline:
- 101mA vs 1a00A: 4 DCCP blocks, all match
- 101mA vs 1binA: 2 DCCP blocks, all match
- 1a87A vs 1allA: 4 DCCP blocks, all match

All 10 DCCP entries matched with **0.0 score difference** vs Fortran output.

### PyO3 Binding Tests

6/6 tests passing: `read_dat`, `ca_coords`, `import_pdb`, `ProteinStore`,
`compare_pair`, `run_wolf_path`.

## 3. Performance

### Benchmarks (leda, Xeon, release + LTO + target-cpu=native + PGO + bit-packed)

| Benchmark | Rust | Fortran | Ratio |
|-----------|------|---------|-------|
| Single pair (101mA vs 1a00A) | 0.631s | 0.39s | 1.6x |
| 20 directed pairs (5 structures) | 19.0s | 7.76s | 2.45x |
| WOLF-only (20 pairs) | 0.30s | — | — |

### Where Time Goes

PARSI (branch-and-bound) accounts for ~96% of pipeline runtime. WOLF+DP+DALICON
is ~0.3s for 20 pairs vs ~19s total.

### Optimization History

| Change | 20-pair time | Cumulative speedup |
|--------|-------------|-------------------|
| Baseline (release defaults) | 25.37s | 1.0x |
| + LTO, HashMap→Vec, stack scoring table | 22.72s | 1.12x |
| + SearchState flat arrays | 22.78s | 1.11x |
| + target-cpu=native, #[inline(always)] | 20.73s | 1.22x |
| + PGO | 20.0s | 1.27x |
| + Bit-packed SearchState + PGO | 19.0s | 1.34x |

Remaining gap vs Fortran (2.45x) is in PARSI scoring inner loops. Further gains
require algorithmic changes (pruning, SIMD), not data structure optimization.

## 4. How DPAM Currently Uses DALI

### Step 7: Iterative DALI (`dpam/steps/step07_iterative_dali.py`)

For each protein, Step 7:
1. Reads ECOD domain candidates from Step 6 (10-100+ domains per protein)
2. Launches `multiprocessing.Pool` workers (one per domain)
3. Each worker runs an **iterative DALI loop**:
   - Calls `dali.pl --pdbfile1 <query> --pdbfile2 <template> --outfmt summary,alignments,transrot`
   - Parses `mol*.txt` output files for z-score, alignment ranges, rotation matrix, translation vector
   - If alignment < 20 residues, stops
   - Removes aligned residues from query PDB (writes new PDB with remaining atoms)
   - Repeats until < 20 residues remain
4. Concatenates all per-domain hit files into `{prefix}_iterativdDali_hits`

### DALI Wrapper (`dpam/tools/dali.py`)

```python
class DALI(ExternalTool):
    def align(self, pdb1, pdb2, output_dir, ...) -> (z_score, alignments, rotation_rows, translation_vals):
        # Spawns subprocess: dali.pl --pdbfile1 ... --pdbfile2 ... --outfmt summary,alignments,transrot
        # Parses mol*.txt output
        # Returns: (float, [(q_resid, t_resid), ...], [rot_row_str, ...], [trans_val_str, ...])
```

Returns a 4-tuple:
- `z_score`: float (or None)
- `alignments`: list of `(query_resid, template_resid)` pairs (1-based)
- `rotation_rows`: list of 3 tab-separated strings (`"r11\tr12\tr13"`)
- `translation_vals`: list of 3 string values (`["t1", "t2", "t3"]`)

### Step 8: Analyze DALI (`dpam/steps/step08_analyze_dali.py`)

Reads the concatenated hits file, computes weighted Q-scores, percentile rankings,
and outputs `{prefix}_good_hits` with 15 columns including rotation/translation.

### Output Format (Step 7 → Step 8)

```
>{edomain}_{iteration}\t{zscore}\t{n_match}\t{qlen}\t{slen}
rotation\tval1\tval2\tval3
rotation\tval1\tval2\tval3
rotation\tval1\tval2\tval3
translation\tval1\tval2\tval3
{query_resid}\t{template_resid}
{query_resid}\t{template_resid}
...
```

### Per-Invocation Overhead in Current Setup

Each `dali.pl` call:
1. **Process creation**: fork/exec Perl interpreter + Fortran binary
2. **PDB → DAT import**: Parse PDB, compute DSSP, build domain tree, write `.dat` files
3. **Pipeline execution**: WOLF → DP → DALICON (+ PARSI + FILTER95 for full run)
4. **Output**: Write `mol*.txt`, log files
5. **Parse**: Python reads and parses text output files

For a ~0.39s comparison, process/IO overhead is a significant fraction.
A typical protein triggers 50-500+ of these invocations across all candidate domains.

### Scale

| Scenario | DALI invocations | Est. time (current) |
|----------|-----------------|---------------------|
| 1 protein, 50 domains, 2 iterations avg | ~100 | ~40s (parallel) |
| 100 proteins | ~10,000 | ~1 hour (8 workers) |
| 1,000 proteins | ~100,000 | ~10 hours |

## 5. What the Rust Library Provides vs What DPAM Needs

### Currently Available

| Feature | Status | Notes |
|---------|--------|-------|
| Full pipeline (WOLF+PARSI both paths) | Yes | `compare_pair()` |
| Individual paths | Yes | `run_wolf_path()`, `run_parsi_path()` |
| PDB/CIF import | Yes | `import_pdb()` — standalone, no external tools |
| .dat file I/O | Yes | `read_dat()` |
| In-memory protein caching | Yes | `ProteinStore` with Arc/RwLock |
| Z-score | Yes | `DccpEntry.zscore` |
| Alignment blocks (ranges) | Yes | `DccpEntry.blocks` → `[AlignmentBlock(l1,r1,l2,r2)]` |
| Score | Yes | `DccpEntry.score` |
| RMSD | Yes | `DccpEntry.rmsd` |

### Missing for DPAM Integration

| Feature | Status | Why DPAM Needs It |
|---------|--------|-------------------|
| **Rotation matrix** | NOT in DccpEntry | Step 7/8 output includes 3x3 rotation for each hit |
| **Translation vector** | NOT in DccpEntry | Step 7/8 output includes 3-component translation |
| **Residue-level alignment pairs** | Partial | DccpEntry has block ranges; Step 7 writes individual `(q_resid, t_resid)` pairs. Blocks can be expanded but this is a convenience gap. |
| **PDB subsetting** | NOT in library | Step 7 removes aligned residues and writes a new PDB for the next iteration. Currently done in Python by filtering ATOM lines. |
| **Single-direction comparison** | NOT exposed | `compare_pair` runs both directions + both paths. Step 7 only needs query→template, WOLF path only. |
| **Import from PDB with residue ID mapping** | Partial | `import_pdb` generates internal 1-based numbering. Step 7 needs to map between PDB residue IDs and DALI's internal numbering. |

### Interface Gap Analysis

The fundamental mismatch: `compare_pair()` runs the **full** pipeline (both
directions, WOLF + PARSI paths, self-comparison filtering) and returns high-level
`DccpEntry` results. DPAM Step 7 needs a **single-direction, single-path** alignment
with **structural transformation** (rotation/translation) output and **residue-level**
alignment details.

The Fortran `dali.pl` is invoked as a pairwise aligner, not an all-vs-all tool.
Our Rust API was designed for the latter use case.

## 6. Build and Runtime Dependencies

### Build-time
- Rust 1.93.1 (stable)
- maturin 1.12
- System: Linux x86_64

### Rust Crates
- `ndarray` 0.16 — matrices
- `pdbtbx` 0.12 — PDB/CIF parsing
- `flate2` 1 — gzip support
- `pyo3` 0.23 — Python bindings
- `numpy` 0.23 — numpy array interop
- `approx` 0.5 — float comparison (tests only)

### Runtime
- Python 3.11 (Anaconda on leda; system Python 3.10 also works)
- `numpy` (Python package)
- `LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libgomp.so.1` on leda for numpy

### Compiler Flags (release profile)
- `opt-level = 3`
- `lto = "fat"`
- `codegen-units = 1`
- `target-cpu = native` (via `.cargo/config.toml`)
- PGO: optional, applied via `RUSTFLAGS` (not in config)

## 7. File Inventory

```
dali_rust/
├── Cargo.toml                    # Workspace: dali-core + dali-python
├── .cargo/config.toml            # target-cpu=native
│
├── dali-core/
│   ├── Cargo.toml                # Library crate
│   ├── src/                      # 9,015 lines, 32 files
│   │   ├── lib.rs                # Module declarations
│   │   ├── store.rs              # ProteinStore (Arc/RwLock caching)
│   │   ├── pipeline.rs           # E2E orchestrator (compare_pair, run_*_path)
│   │   ├── types/                # Protein, Segment, AlignmentBlock, DccpEntry
│   │   ├── numerics/             # nint, Kabsch/u3b, fitz, NW, scoring
│   │   ├── io/                   # .dat parser, PDB reader, DSSP, secstr, domain, import
│   │   ├── wolf/                 # SSE spatial hashing (4 files)
│   │   ├── dp/                   # Z-score computation (1 file)
│   │   ├── dalicon/              # GA refinement (1 file)
│   │   ├── parsi/                # Branch-and-bound (6 files)
│   │   └── filter95/             # Redundancy filtering + FITZ (1 file)
│   │
│   └── tests/                    # 2,075 lines, 10 files, 60 tests
│       ├── test_wolf.rs          # 4 tests (5-struct + 18-struct groups)
│       ├── test_dp.rs            # 2 tests
│       ├── test_dalicon.rs       # 4 tests
│       ├── test_parsi.rs         # 4 tests
│       ├── test_filter95.rs      # 4 tests
│       ├── test_e2e.rs           # 1 test (3 reference pairs, 10 entries)
│       ├── test_import.rs        # 4 tests (PDB reader + import pipeline)
│       ├── test_dat_reader.rs    # 10 tests
│       ├── test_store.rs         # 4 tests
│       └── test_scoring.rs       # 2 tests
│
└── dali-python/
    ├── Cargo.toml                # cdylib crate
    ├── pyproject.toml            # maturin config
    ├── src/lib.rs                # 213 lines — PyO3 wrappers
    ├── test_bindings.py          # 6 Python tests
    └── bench.py                  # Performance benchmark script
```

## 8. Known Limitations

1. **No rotation/translation in output.** The Kabsch superposition (u3b) is computed
   internally during DALICON/FITZ but the rotation matrix and translation vector are
   not propagated to DccpEntry. They are discarded after scoring.

2. **No single-direction API.** `compare_pair` runs both cd1→cd2 and cd2→cd1. For
   DPAM's iterative alignment, only one direction is needed per call.

3. **No PDB subsetting.** The iterative loop in Step 7 removes aligned residues from
   the query PDB and re-runs DALI. The current Rust API has no mechanism for this —
   it operates on whole proteins loaded from .dat files or PDB files.

4. **ProteinStore is .dat-file-oriented.** It caches proteins keyed by code, loaded
   from `{dat_dir}/{code}.dat`. DPAM Step 7 works with PDB files that change each
   iteration (residues removed). The store would need to support ad-hoc Protein
   objects, not just file-backed ones.

5. **Residue numbering.** DALI uses internal 1-based sequential numbering (1..nres).
   PDB files have their own residue numbering (which may have gaps, insertion codes).
   Step 7 maps between these. The Rust `import_pdb` creates sequential numbering but
   doesn't expose the PDB→internal mapping.

6. **Thread safety for multiprocessing.** `ProteinStore` uses `Arc<RwLock>` which is
   thread-safe within a single process. Python's `multiprocessing.Pool` (used by
   Step 7) forks separate processes — each would need its own `ProteinStore`. However,
   template proteins could be shared via a read-only store loaded before forking.

7. **PARSI match rate.** The Rust PARSI matches Fortran 71.6% (expanded corpus) vs
   Python's 78%. This is accumulated float32 search divergence, not a correctness
   issue. End-to-end results are unaffected (0.0 score difference on reference pairs).

## 9. What Integration Requires

To replace `dali.pl` subprocess calls in DPAM Step 7, we need:

1. **Add rotation/translation to DccpEntry** — propagate u3b results through the pipeline
2. **Expose a pairwise alignment function** — single direction, returns z-score + alignments + transform
3. **Support ad-hoc Protein objects** — allow Step 7 to create/modify Protein structs (e.g., subset residues) without writing .dat files
4. **Residue ID mapping** — expose PDB residue numbering alongside internal numbering
5. **Expand alignment blocks to residue pairs** — or let Step 7/8 do it (blocks are ranges, not individual pairs)

The largest architectural change is (3): DPAM's iterative loop mutates the query
each iteration. Options:
- **A.** Re-import from PDB each iteration (Step 7 already writes a new PDB file)
- **B.** Provide a Protein subset operation in Rust
- **C.** Skip iteration entirely — run full DALI once, return all alignment blocks, let Step 8 handle decomposition

Option A is simplest and matches current Step 7 behavior. `import_pdb` takes ~5ms
(vs ~400ms for a comparison), so the overhead is negligible.
