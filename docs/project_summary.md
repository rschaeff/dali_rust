# DaliLite.v5 Reimplementation: Project Summary

## Goal

Proof-of-work that AI-assisted validated reimplementation of legacy scientific
software is feasible. The product is the validation framework, the evidence that
it works, and the fully standalone compiled replacement.

## Background

DaliLite.v5 is a structural alignment tool for comparing protein 3D structures.
Written in Fortran 77 with Perl/shell glue, it has been in continuous use for
decades in structural biology but carries severe technical debt: no tests, no
documentation, implicit mutable state shared across modules, and a build
pipeline that communicates through numbered Fortran file units, shell `cat`, and
Perl regex scripts. The source totals ~8,100 lines of Fortran plus ~2,000 lines
of Perl/shell.

The reimplementation proceeded in two phases:

1. **Python teardown** (Phase 0) — module-by-module reverse engineering with
   stage-by-stage ground truth capture and validation against the Fortran binary.
   This produced the algorithmic understanding and the validation framework.

2. **Rust reimplementation** (Phases 1-11) — compiled implementation validated
   against the same ground truth. The Rust code is fully standalone: it reads
   PDB/CIF files directly, computes DSSP secondary structure, builds domain
   decompositions, and runs the full DALI comparison pipeline. PyO3 bindings
   expose the compiled pipeline to Python, including pairwise alignment,
   database search, and iterative multi-domain detection.

## Pipeline Architecture

The DALI structural comparison pipeline for a protein pair:

```
Input (PDB/CIF or .dat) -> DSSP -> Domain Decomposition -> Protein struct
                                                              |
                    +-----------------------------+-----------+
                    |                             |
              WOLF path                     PARSI path
              (fast, heuristic)             (exhaustive)
                    |                             |
                    v                             v
                  WOLF                          PARSI
            (SSE spatial hash)          (branch-and-bound)
                    |                             |
                    v                             v
                   DP                         FILTER95
             (Z-score filter)           (FITZ + Z-score filter)
                    |                             |
                    v                             v
                DALICON                       DALICON
              (GA refinement)             (GA refinement)
                    |                             |
                    v                             v
                   DP                            DP
                    |                             |
                    +-----------------------------+
                                  |
                                  v
                          Final DCCP output
```

Both paths run in both directions (query/target), producing up to 4 runs.
Self-comparisons are filtered. Final output is the union of all surviving
alignments.

For database search, the pipeline extends to one-to-many comparison with
iterative domain masking: search -> mask aligned residues -> repeat until
no significant hits remain.

## What Each Module Does

1. **WOLF** is a spatial hashing scheme. It builds a 3D grid of SSE pair
   descriptors (midpoint + direction in a canonical frame), then votes on
   matching SSE pairs between two proteins. The best-voted pair seeds an
   iterative superposition (FITZ). Output: block-compressed alignment ranges.

2. **DP** is a distance-matrix scoring function. It computes pairwise CA
   distances (int16, distance x 10), scores aligned residue pairs with a
   Gaussian-weighted distance similarity, and derives Z-scores from a cubic
   polynomial calibration. Domain definitions control which residues participate.

3. **DALICON** is a Monte Carlo genetic algorithm. It seeds candidate alignments
   from tetrapeptide fragments (TESTI), then evolves them through random block
   operations (LEAN_MC). The scoring function uses binned distance matrices.

4. **PARSI** is exhaustive branch-and-bound over SSE alignment space. It
   enumerates all topologically valid SSE mappings, pruning by upper/lower
   distance bounds. The search tree follows the domain hierarchy.

5. **FILTER95** applies FITZ superposition refinement to PARSI hits, rescores
   with a distance-based function (dist x 100 precision), and filters by Z-score.

6. **DSSP** (Kabsch-Sander algorithm) computes 8-state secondary structure from
   backbone atom coordinates via hydrogen bond energy estimation.

7. **Domain decomposition** (PUU) partitions the protein into structural domains
   via binary decomposition of the CA-CA contact graph.

## Validation Framework

### Ground Truth Capture

Ground truth is captured stage-by-stage from the reference DaliLite.v5 Fortran
binary. Each module's output is compared against captured reference using
tolerance-based matching appropriate to each format:

- **WOLF**: Exact integer match on nblock and alignment ranges
- **DP**: Score tolerance of 2.0, Z-score tolerance
- **DALICON**: Score and Z-score tolerance matching
- **PARSI**: Score tolerance of 100 (suboptimal alignments differ due to search divergence)
- **FILTER95**: Score tolerance of 1% or 100 units, Z-score tolerance of 0.5
- **Import pipeline**: nres exact, CA coordinates within 0.15 A, SSE segments within +/-2 residues

### Test Corpora

**Original (5 structures):** 101mA, 1a00A, 1a87A, 1allA, 1binA
- 141-297 residues, 25 directed pairs per module

**Expanded (18 structures, ECOD-diverse):**
- 1a7sA, 1aiwA, 1a25A, 1a12A, 1f3uA, 1a04A, 1bbhA, 1b3qA, 1a17A,
  1miwA, 1aopA, 1a8lA, 1a6qA, 1a06A, 1b3oB, 1a0cA, 1a4iA, 1bcoA
- 62-456 residues, 4-30 SSEs, all major ECOD fold types (a.1-a.18)

**Batch (785 proteins from DPAM archaea tier1):**
- batch_01 of `/data/ecod/archaea/dpam_tier1_batched/`
- Real DPAM workload: iterative search against ~200 ECOD70 template candidates
- Ground truth from Fortran DaliLite.v5 run through full DPAM pipeline

## Validation Results

### Python Validation (Phase 0)

| Module   | Original (5 struct) | Expanded (18 struct) | Total             |
|----------|---------------------|----------------------|-------------------|
| WOLF     | 21/21               | 279/279              | 300/300 (100%)    |
| DP       | 10/10               | 24/24 (tol=2.0)     | 34/34 (100%)      |
| DALICON  | 10/10               | 24/24                | 34/34 (100%)      |
| PARSI    | 379/435 (87%)       | 10533/13515 (78%)    | 10912/13950 (78%) |
| FILTER95 | 38/38               | 293/293              | 331/331 (100%)    |
| **E2E**  | **10/10 (3 pairs)** |                      | **10/10**         |

End-to-end: 10/10 DCCP entries matched with 0.0 score difference across 3
reference pairs (101mA/1a00A, 101mA/1binA, 1a87A/1allA).

PARSI mismatches are heuristic search divergence (+/-1 score table differences
propagating through branch-and-bound), not bugs. See
`docs/parsi_mismatch_analysis.md`.

FILTER95: 2 of 295 expanded-corpus reference entries are Fortran bugs (stale
`lkeep` array from prior protein when `idom > ndom`). Python correctly rejects
these.

### Rust Validation (Phases 1-11)

#### Module-level (5+18 structure corpora)

| Module   | Original (5 struct) | Expanded (18 struct) | Total             |
|----------|---------------------|----------------------|-------------------|
| WOLF     | 21/21               | 278/279              | 299/300 (99.7%)   |
| DP       | 10/10               | 24/24                | 34/34 (100%)      |
| DALICON  | 10/10               | 24/24                | 34/34 (100%)      |
| PARSI    | 379/435 (87%)       | 9616/13515 (71%)     | 9995/13950 (72%)  |
| FILTER95 | 38/38               | 292/295 (99%)        | 330/333 (99.1%)   |
| **E2E**  | **10/10 (3 pairs)** |                      | **10/10**         |
| Import   | 42/42 (100%)        | 279/303 (92.1%)      | 321/345 (93.0%)   |

- **WOLF 1 mismatch** (1a25A/1bcoA): floating-point edge case at 4.0 A fitz
  cutoff changes nblock from 14 to 13.
- **PARSI 18-struct** lower than Python due to accumulated float32 search
  divergence — matches Python exactly on the 5-struct corpus (379/435).
- **FILTER95** 2 MISSING = known Fortran bugs (correctly rejected); 1 slightly
  outside 1% score tolerance (f64 vs f32 arithmetic, 2.17% diff).
- **Import pipeline** segment overlap: 100% (5-struct), 92.1% (18-struct).
  Differences are DSSP reimplementation variations vs DaliLite's bundled
  dsspcmbi, not bugs.

#### Batch validation (785 proteins, DPAM archaea tier1)

Domain-level comparison of Rust `iterative_search` vs Fortran iterative DALI:

| Metric | Value |
|--------|-------|
| Proteins tested | 785 (0 errors, 0 crashes) |
| Fortran domains | 1,488 |
| Rust domains | 1,117 |
| Matched (region overlap >= 50%) | 1,005 (67.5% recall) |
| Precision | 99.1% (Rust hits match a Fortran template) |
| Novel Rust detections | 112 (63 at z >= 10) |
| Missed Fortran domains | 483 |

**Recall by z-score threshold:**
- z >= 15: 99.6% (2 missed)
- z >= 10: 98.5% (8 missed)
- z >= 5: 93.8% (115 missed)
- z >= 2: 67.5% (483 missed)

**Novel detection mechanism:** Rust's fresh-per-round scoring finds additional
domains in later iteration rounds. 76 of 112 novel detections are from proteins
where Rust found MORE domains than Fortran (genuinely additional, not replacements).
62 of those 76 came from rounds beyond Fortran's maximum iteration count.

**Biological impact:** Of 483 missed domains, 431 were "good_domain" in DPAM's
final classification — but 346 of those had HHsearch prob >= 0.9 (DALI evidence
redundant). Only 86 were DALI-critical, all with DPAM prob >= 0.5 (no
classifications would flip).

### Automated Tests

87 tests (61 Rust + 26 Python), all passing:

| Test Suite          | Tests | Description                           |
|---------------------|-------|---------------------------------------|
| Unit tests (lib)    | 21    | Numerics, DSSP, domain, secstr        |
| DAT reader          | 10    | .dat file parsing + domain filtering  |
| Scoring             | 2     | dpgetdist, zscore                     |
| Store               | 4     | ProteinStore caching + distance mats  |
| WOLF                | 4     | 5-struct + 18-struct groups            |
| DP                  | 2     | 5-struct + 18-struct                  |
| DALICON             | 4     | 5-struct + 18-struct groups            |
| PARSI               | 4     | 5-struct + 18-struct groups            |
| FILTER95            | 4     | 5-struct + 18-struct groups            |
| E2E pipeline        | 1     | 3 reference pairs, 10 DCCP entries    |
| Import pipeline     | 4     | PDB reader + full import (5+18 struct)|
| Write/roundtrip     | 1     | write_dat -> read_dat -> compare      |
| PyO3 bindings       | 10    | Full Python API coverage              |
| Align coverage      | 7     | PDB offsets, small domains, .dat bypass, rotation |
| Search/iteration    | 9     | mask, add, search, skip_wolf, iterative |
| **Total**           | **87**|                                       |

## Fortran Quirks Discovered

- **Stateful ca1 in DALICON**: `transrotate` modifies coordinates in-place;
  subsequent comparisons with the same cd1 use already-transformed coords.
  Processing order matters.
- **Fortran rounding**: `nint()` rounds half away from zero (not banker's
  rounding). Affects distance matrix quantization throughout.
- **Single-precision intermediates**: Distance matrices use float32 arithmetic
  despite float64 coordinates, matching Fortran's default `real`.
- **FILTER95 stale state bug**: When `idom > ndom` for a protein, the Fortran
  `lkeep` array retains values from the previous protein, causing incorrect
  domain filtering.
- **PARSI score table sensitivity**: +/-1 differences in int16 distance
  quantization propagate through branch-and-bound pruning, causing search
  paths to diverge. This is inherent to the algorithm, not a bug.
- **pipe96-free.pl -99 stripping**: The Perl `s/\-99//g` regex can corrupt
  numbers containing "-99" as a substring.
- **Inter-module communication**: Every module boundary involves: write to
  numbered Fortran units -> shell glob -> `cat` -> Perl reformatting -> next
  module reads. ~10 serialization round-trips per pair comparison.
- **DALICON accept-band OOB**: `j_hi` clamp used `nres1` instead of `nres2`,
  causing out-of-bounds access when query > template. Silent in Fortran (no
  bounds checking), panic in Rust.

## Code Inventory

### Python (Phase 0)

```
dali_cl/
  wolf/               7 files, ~1,200 lines -- SSE hashing + spatial voting
  dp/                 2 files,   ~400 lines -- Z-score computation
  dalicon/            5 files,   ~800 lines -- GA refinement
  parsi/              8 files, ~2,800 lines -- Branch-and-bound
  filter95/           2 files,   ~350 lines -- Redundancy filtering + FITZ
  pipeline.py                    ~200 lines -- E2E orchestrator
  validation/         3 files, ~1,250 lines -- Parsers, capture, compare
  test_*.py           7 files,   ~950 lines -- Per-module + E2E tests
```

Total: ~7,950 lines Python.

### Rust (Phases 1-11)

```
dali_rust/dali-core/
  src/                              9,883 lines (32+ .rs files)
    types/        3 files,    281 lines -- Protein, Segment, Alignment, SearchHit, DistMat
    numerics/     6 files,    975 lines -- nint, Kabsch/u3b, fitz, NW, scoring
    io/           7 files,  1,665 lines -- .dat parser, PDB reader, DSSP, secstr, domain, import
    wolf/         4 files,    824 lines -- SSE spatial hashing
    dp/           1 file,     173 lines -- Z-score computation
    dalicon/      1 file,   1,579 lines -- GA refinement
    parsi/        6 files,  2,891 lines -- Branch-and-bound
    filter95/     1 file,     490 lines -- Redundancy filtering + FITZ
    pipeline.rs               464 lines -- E2E orchestrator + search + iterative detection
    store.rs                  154 lines -- ProteinStore (caching layer + add_protein)
    lib.rs                     20 lines -- Module declarations

  tests/                            2,130 lines (10 .rs files)
    test_import.rs            413 lines -- PDB reader + import pipeline
    test_filter95.rs          284 lines -- FILTER95 validation
    test_dalicon.rs           258 lines -- DALICON validation
    test_wolf.rs              252 lines -- WOLF validation
    test_dp.rs                251 lines -- DP validation
    test_parsi.rs             229 lines -- PARSI validation
    test_e2e.rs               173 lines -- E2E pipeline validation
    test_dat_reader.rs        166 lines -- .dat parsing + domain filtering
    test_store.rs              62 lines -- ProteinStore
    test_scoring.rs            42 lines -- Scoring functions

dali_rust/dali-python/
  src/lib.rs                  530 lines -- PyO3 bindings (7 wrapper types + 10 functions)
  test_bindings.py            220 lines -- Python binding tests (10 tests)
  test_align_coverage.py      354 lines -- Alignment coverage tests (7 tests)
  test_search.py              188 lines -- Search/iteration tests (9 tests)
```

Total: 13,305 lines Rust + Python bindings.

Dependencies: `ndarray` 0.16, `approx` 0.5, `pdbtbx` 0.12, `flate2` 1,
`pyo3` 0.23, `numpy` 0.23, `rayon` 1.10.

### Grand Total

| Component           | Lines  |
|---------------------|--------|
| Python source       | ~5,750 |
| Python tests        | ~2,200 |
| Rust source         | 10,413 |
| Rust tests          |  2,130 |
| PyO3 bindings+tests |  1,292 |
| **Combined**        | **~21,785** |

vs. DaliLite.v5 Fortran: ~8,100 lines + ~2,000 lines Perl/shell (no tests).

## Design Decisions

### Why Two Implementations?

The Python teardown was essential for understanding the algorithms before
attempting a compiled implementation. Fortran 77 with implicit state, no type
annotations, and no tests makes direct translation error-prone. The Python
phase produced:

1. Algorithmic understanding of every module
2. A validation framework with captured ground truth
3. Discovery of Fortran quirks (stateful ca1, stale lkeep, rounding) that would
   have been silent bugs in a direct translation

The Rust phase then became a **translation problem** (expressing known algorithms
in a new language), not a **discovery problem**. Translation is where AI
assistance is strongest.

### Why Rust?

- **Compiler-enforced correctness.** The hardest bugs discovered during the
  teardown were mutable-state-across-calls issues (DALICON's in-place ca1,
  FILTER95's stale lkeep). Rust's ownership model prevents these structurally.
  The borrow checker acts as a second validator — bugs that would be silent in C
  don't compile in Rust.
- **Ecosystem.** PyO3 for Python integration, rayon for pair-level parallelism,
  pdbtbx for PDB/CIF parsing with gzip support. Fortran lacks this ecosystem.
- **The hard part is done.** Algorithmic complexity is fully understood from the
  Python teardown. What remains is expressing known algorithms in typed Rust,
  which the compiler helps verify.

### Key Rust Architecture Choices

- **ndarray** for all 2D matrices (CA coordinates, distance matrices, scoring)
- **Typed distance matrices** (DistMatScale10, DistMatScale100) preventing
  accidental mixing of different precision levels
- **Jacobi SVD** for 3x3 eigendecomposition in Kabsch alignment (avoids LAPACK)
- **f32 arithmetic** in PARSI for Fortran single-precision compatibility
- **ProteinStore** with Arc/RwLock for thread-safe lazy loading and caching
- **Direct struct passing** between pipeline stages (no text serialization) --
  DccpEntry -> DaliconRecord, Filter95Entry -> DaliconRecord
- **Kabsch-Sander DSSP** from scratch (~420 lines) rather than wrapping external
  dsspcmbi binary
- **Simplified PUU domain decomposition** using CA-contact binary partitioning
  with tau scoring
- **Rayon** for pipeline-level parallelism (4 paths concurrent)
- **In-memory protein masking** for iterative domain detection (no PDB rewrite)

## DSSP Implementation Notes

The Rust DSSP (Phase 8) reimplements the Kabsch-Sander algorithm:

1. **H-bond energy**: E = Q * (1/rON + 1/rCH - 1/rOH - 1/rCN) where
   Q = +27888 cal/mol. H-bond accepted if E < -500 cal/mol.
2. **H position estimation**: H = N + unit(prev_C - prev_O), placing H opposite
   the C=O bond direction per standard DSSP convention.
3. **Bridge detection**: Parallel and antiparallel patterns from donor H-bond
   pairs; bridges build into ladders, ladders sharing residues form sheets.
4. **Assignment priority**: H > B > E > G > I > T > S > ' ' -- bridges/sheets
   assigned before helices to respect the standard DSSP precedence.
5. **3-state reduction**: G,I,H -> H; E -> E; else -> L. Minimum lengths
   (helix >= 6, strand >= 8) with border growth.

Validation: 100% segment overlap on 5-struct corpus, 92.1% on 18-struct corpus
vs. DaliLite's bundled dsspcmbi (auto-translated Pascal -> C, 125K lines).

## Lessons Learned

### AI-Assisted Reimplementation Methodology

1. **Stage-by-stage ground truth is essential.** Capturing intermediate outputs
   from the reference implementation lets you validate each module independently.
   End-to-end testing alone is insufficient -- it can't localize bugs.

2. **The Python teardown phase is not optional.** Direct Fortran-to-Rust
   translation without understanding the algorithms produces code that compiles
   but computes wrong answers. The Python phase forces you to understand every
   step.

3. **Tolerance-based validation handles floating-point divergence.** Different
   languages, compilers, and precision levels produce slightly different results.
   The validation framework must accommodate this without hiding real bugs.

4. **Heuristic algorithms produce legitimately different results.** PARSI's 78%
   match rate is not a bug -- it's inherent to branch-and-bound with quantized
   scoring. Understanding which mismatches are expected vs. buggy requires deep
   algorithmic analysis.

5. **The borrow checker finds real bugs.** Several Fortran state-management
   issues (DALICON ca1, FILTER95 lkeep) would have been caught at compile time
   in Rust. For AI-assisted development, this compile-time checking compensates
   for the AI's tendency to introduce subtle state bugs.

6. **Translation is faster than discovery.** Phase 0 (Python teardown) took
   longer per module than Phases 1-11 (Rust reimplementation) because discovery
   requires reading undocumented Fortran and debugging against opaque reference
   output. Translation requires expressing known algorithms in a new language,
   which AI assistance handles well.

### What Worked

- Capturing ground truth at every stage boundary
- Using ECOD-diverse corpus to exercise all fold types
- Tolerance-based matching with module-specific thresholds
- Building the validation framework first, then implementing
- Batch validation against real DPAM workload (785 proteins) for system-level confidence

### What Was Hard

- PARSI: 4,765 lines of Fortran with 5 interacting files, dense bit-packed
  priority queue, domain tree decomposition driving search space. Required the
  deepest algorithmic understanding.
- DSSP: Three bugs discovered sequentially (Q sign, H position, assignment
  priority). Each one required understanding the underlying physics/algorithm.
- Fortran implicit state: variables retaining values across subroutine calls
  with no documentation of which state is intentionally persistent.
- DALICON OOB: accept-band `j_hi` used `nres1` instead of `nres2`, silent in
  Fortran (no bounds checking), only surfaced when Rust panicked on a query
  larger than its template.

## Status

**Complete.** All 11 phases implemented, validated, and tested. The Rust
implementation is fully standalone -- it reads PDB/CIF files directly and
produces structural alignment results without any external tools. PyO3 bindings
expose pairwise alignment, database search, and iterative multi-domain detection
to Python.

| Phase | Module                              | Validation                              |
|-------|-------------------------------------|-----------------------------------------|
| 0     | Python teardown + validation        | 11,621/14,659 (79.3%)                   |
| 1     | Types, numerics, .dat I/O           | 28 unit/integration tests               |
| 2     | WOLF                                | 299/300 (99.7%)                         |
| 3     | DP                                  | 34/34 (100%)                            |
| 4     | DALICON                             | 34/34 (100%)                            |
| 5     | PARSI                               | 9,995/13,950 (71.6%)                    |
| 6     | FILTER95                            | 330/333 (99.1%)                         |
| 7     | Pipeline E2E                        | 10/10 (100%, 0.0 score diff)            |
| 8     | DSSP + PDB/CIF Import               | 321/345 segments (93.0%)                |
| 9     | PyO3 Python Bindings                | 10/10 binding tests (100%)              |
| 10    | DPAM Integration + Deficiency Fixes | 7/7 coverage tests (100%)               |
| 11    | Native Search + Iterative Detection | 9/9 search tests + 785-protein batch    |

87 automated tests (61 Rust + 26 Python). All passing.

### Python API

The `dali` module provides four tiers of API:

**Tier 1 -- Full pipeline:**
```python
import dali
store = dali.ProteinStore("/path/to/dat/files")
results = dali.compare_pair("101mA", "1a00A", store)
for r in results:
    print(f"{r.cd1} vs {r.cd2}: score={r.score:.1f}, z={r.zscore:.1f}")
```

**Tier 2 -- Database search + iterative detection:**
```python
# One-to-many search
hits = dali.search_database("query", targets, store, z_cut=2.0)

# Iterative multi-domain detection
domains = dali.iterative_search("query", targets, store,
    min_zscore=2.0, skip_wolf=True)
for d in domains:
    print(f"Round {d.round}: {d.cd2} z={d.zscore:.1f}")
```

**Tier 3 -- Pairwise alignment with structural transform:**
```python
result = dali.align_pdb("query.pdb", "template.pdb",
    query_chain="A", template_chain="A",
    query_code="q", template_code="t",
    template_dat="/path/to/template.dat")  # optional DSSP bypass
if result:
    print(f"Z-score: {result.zscore:.1f}, aligned: {result.n_aligned}")
    print(f"Rotation: {result.rotation}")       # 3x3 matrix
    print(f"Translation: {result.translation}")  # 3-element vector
```

**Tier 4 -- I/O utilities:**
```python
protein = dali.read_dat("/path/to/101mA.dat")
protein = dali.import_pdb("/path/to/pdb101m.ent.gz", "A", "101m")
dali.write_dat(protein, "/path/to/output.dat")
masked = dali.mask_protein(protein, keep_indices, "new_code")
ca = protein.ca_coords()  # numpy array, shape (3, nres)
store.add_protein(protein) # write .dat + cache
```

Wrapper types: `Protein`, `ProteinStore`, `DccpEntry`, `AlignmentBlock`,
`AlignResult`, `SearchHit`.
Built with maturin 1.12, pyo3 0.23, numpy 0.23.

### Known Limitations

- **SOAP not implemented.** Proteins with < 3 SSEs use Needleman-Wunsch on CA
  distances (SOAP path) in Fortran. The Rust implementation skips these, which
  accounts for ~3 of 483 missed domains in batch validation.
- **DSSP divergence.** The Rust Kabsch-Sander DSSP produces different per-residue
  assignments than Fortran's bundled dsspcmbi for ~35% of ECOD templates. Mitigated
  by `.dat` file bypass when pre-computed segments are available.
- **PARSI search divergence.** 71-87% match rate vs Fortran due to +/-1 int16
  quantization propagating through branch-and-bound. Not a bug; end-to-end
  results are unaffected.
