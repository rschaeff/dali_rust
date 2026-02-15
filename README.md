# DaliLite.v5 Reimplementation

A validated reimplementation of the [DaliLite.v5](http://ekhidna2.biocenter.helsinki.fi/dali/) structural alignment tool in Python and Rust, with PyO3 Python bindings. This project demonstrates that AI-assisted reimplementation of legacy scientific software is feasible when paired with rigorous stage-by-stage validation.

DaliLite.v5 compares protein 3D structures by aligning their distance matrices. Written in ~8,100 lines of Fortran 77 with ~2,000 lines of Perl/shell glue, it has been in continuous use for decades in structural biology but carries severe technical debt: no tests, no documentation, implicit mutable state shared across modules, and a build pipeline that communicates through numbered Fortran file units, shell `cat`, and Perl regex scripts.

This reimplementation provides a fully standalone compiled replacement that reads PDB/CIF files directly, computes DSSP secondary structure, builds domain decompositions, and runs the full DALI comparison pipeline with no external dependencies.

## Pipeline Architecture

```
Input (PDB/CIF or .dat) --> DSSP --> Domain Decomposition --> Protein
                                                                |
                    +---------------------------+---------------+
                    |                           |
              WOLF path                    PARSI path
              (fast, heuristic)            (exhaustive)
                    |                           |
                  WOLF                        PARSI
            (SSE spatial hash)          (branch-and-bound)
                    |                           |
                   DP                       FILTER95
             (Z-score filter)         (FITZ + Z-score filter)
                    |                           |
                DALICON                     DALICON
              (GA refinement)           (GA refinement)
                    |                           |
                   DP                          DP
                    |                           |
                    +---------------------------+
                                |
                        Final DCCP output
```

Both paths run in both directions (query <-> target). Self-comparisons are filtered. Final output is the union of all surviving alignments.

## Validation Results

Ground truth is captured stage-by-stage from the reference Fortran binary and compared using tolerance-based matching appropriate to each module.

### Test Corpora

- **5-structure corpus**: 101mA, 1a00A, 1a87A, 1allA, 1binA (141-297 residues)
- **18-structure corpus** (ECOD-diverse): 62-456 residues, 4-30 SSEs, all major fold types

### Rust vs Fortran Ground Truth

| Module   | 5-struct | 18-struct | Total | Notes |
|----------|----------|-----------|-------|-------|
| WOLF     | 21/21    | 278/279   | 299/300 (99.7%) | 1 float edge case at fitz cutoff |
| DP       | 10/10    | 24/24     | 34/34 (100%) | |
| DALICON  | 10/10    | 24/24     | 34/34 (100%) | |
| PARSI    | 379/435  | 9,616/13,515 | 9,995/13,950 (71.6%) | Search divergence, not bugs |
| FILTER95 | 38/38    | 292/295   | 330/333 (99.1%) | 2 Fortran bugs correctly rejected |
| Import   | 42/42    | 279/303   | 321/345 (93.0%) | DSSP reimplementation variations |
| **E2E**  | **10/10** |          | **10/10 (100%)** | **0.0 score difference** |

PARSI mismatches are inherent to branch-and-bound with quantized scoring: +/-1 differences in int16 distance quantization propagate through pruning, causing search paths to diverge. End-to-end results are unaffected.

### Automated Tests

61 Rust tests + 9 Python binding tests, all passing:

| Test Suite | Tests | Description |
|-----------|-------|-------------|
| Unit tests | 21 | Numerics, DSSP, domain, secstr |
| DAT reader | 10 | .dat file parsing + domain filtering |
| Scoring | 2 | Distance matrices, Z-scores |
| Store | 4 | ProteinStore caching |
| WOLF | 3 | 5-struct + 18-struct validation |
| DP | 2 | 5-struct + 18-struct |
| DALICON | 4 | 5-struct + 18-struct |
| PARSI | 4 | 5-struct + 18-struct |
| FILTER95 | 4 | 5-struct + 18-struct |
| E2E pipeline | 1 | 3 reference pairs, 10 DCCP entries |
| Import | 4 | PDB reader + full import pipeline |
| Write/roundtrip | 1 | write_dat -> read_dat -> compare |
| PyO3 bindings | 9 | Full Python API coverage |

## Performance

Benchmarked on Xeon (release + LTO + target-cpu=native + PGO + bit-packed):

| Benchmark | Rust | Fortran | Ratio |
|-----------|------|---------|-------|
| Single pair (101mA vs 1a00A) | 0.63s | 0.39s | 1.6x |
| 20 directed pairs (5 structures) | 19.0s | 7.76s | 2.45x |
| WOLF-only (20 pairs) | 0.30s | -- | -- |

PARSI (branch-and-bound) accounts for ~96% of pipeline runtime. The remaining gap vs Fortran is in PARSI scoring inner loops where gfortran benefits from column-major array layout and aggressive scalar optimization.

## Python API

### Installation

```bash
cd dali_rust/dali-python
pip install maturin
maturin develop --release
```

### Usage

```python
import dali

# Full pipeline comparison
store = dali.ProteinStore("/path/to/dat/files")
results = dali.compare_pair("101mA", "1a00A", store)
for r in results:
    print(f"{r.cd1} vs {r.cd2}: score={r.score:.1f}, z={r.zscore:.1f}, rmsd={r.rmsd:.1f}")

# Individual paths
wolf_results = dali.run_wolf_path("101mA", ["1a00A", "1binA"], store)
parsi_results = dali.run_parsi_path("101mA", ["1a00A"], store)

# PDB-to-PDB alignment (imports, compares, returns rotation/translation)
result = dali.align_pdb("query.pdb", "template.pdb", query_chain="A", template_chain="A")
if result:
    print(f"Z-score: {result.zscore:.1f}, aligned: {result.n_aligned}")
    print(f"Rotation: {result.rotation}")       # 3x3 matrix
    print(f"Translation: {result.translation}")  # 3-element vector
    for q_resid, t_resid in result.alignments:   # PDB residue numbering
        print(f"  {q_resid} <-> {t_resid}")

# I/O utilities
protein = dali.read_dat("/path/to/101mA.dat")
protein = dali.import_pdb("/path/to/pdb101m.ent.gz", "A", "101m")
ca = protein.ca_coords()   # numpy array, shape (3, nres)
protein.write_dat("/path/to/output.dat")
```

### Types

| Type | Properties |
|------|-----------|
| `Protein` | `.code`, `.nres`, `.nseg`, `.na`, `.nb`, `.sequence`, `.resid_map`, `.ca_coords()`, `.write_dat()` |
| `ProteinStore` | `__init__(dat_dir)`, `.get_protein(code)`, `__contains__`, `__len__` |
| `DccpEntry` | `.cd1`, `.cd2`, `.score`, `.zscore`, `.rmsd`, `.nblock`, `.blocks` |
| `AlignmentBlock` | `.l1`, `.r1`, `.l2`, `.r2` (1-based residue ranges) |
| `AlignResult` | `.zscore`, `.score`, `.rmsd`, `.n_aligned`, `.alignments`, `.rotation`, `.translation`, `.blocks` |

## Building from Source

### Prerequisites

- Rust 1.75+ (stable)
- Python 3.10+
- numpy

### Rust library

```bash
cd dali_rust
cargo test --release           # Run all 61 tests
cargo build --release          # Build library only
```

### Python bindings

```bash
cd dali_rust/dali-python
pip install maturin
maturin develop --release      # Build + install in current Python env
python test_bindings.py        # Run 9 binding tests
```

## Project Structure

```
dali_cl/
  wolf/               7 files, ~1,200 lines  -- Python: SSE spatial hashing
  dp/                 2 files,   ~400 lines  -- Python: Z-score computation
  dalicon/            5 files,   ~800 lines  -- Python: GA refinement
  parsi/              8 files, ~2,800 lines  -- Python: branch-and-bound
  filter95/           2 files,   ~350 lines  -- Python: redundancy filtering
  pipeline.py                    ~200 lines  -- Python: E2E orchestrator
  validation/         3 files, ~1,250 lines  -- Parsers, capture, compare framework
  test_*.py           7 files,   ~950 lines  -- Per-module + E2E validation

  dali_rust/
    dali-core/
      src/             32 files, ~9,000 lines -- Rust library
        types/         Protein, Segment, AlignmentBlock, DccpEntry
        numerics/      nint, Kabsch/u3b, fitz, NW, scoring
        io/            .dat parser, PDB reader, DSSP, secstr, domain, import
        wolf/          SSE spatial hashing
        dp/            Z-score computation
        dalicon/       GA refinement
        parsi/         Branch-and-bound
        filter95/      Redundancy filtering + FITZ
        pipeline.rs    E2E orchestrator
        store.rs       ProteinStore (thread-safe caching)
      tests/           10 files, ~2,100 lines -- Validation tests

    dali-python/
      src/lib.rs       ~370 lines  -- PyO3 bindings
      test_bindings.py ~150 lines  -- Python binding tests
      bench.py                     -- Performance benchmarks

  docs/                Analysis and project documentation
```

### Dependencies

| Crate | Version | Purpose |
|-------|---------|---------|
| ndarray | 0.16 | 2D matrices (CA coordinates, distance matrices) |
| pdbtbx | 0.12 | PDB/CIF parsing with gzip support |
| flate2 | 1 | Gzip decompression |
| pyo3 | 0.23 | Python bindings |
| numpy | 0.23 | Numpy array interop |
| approx | 0.5 | Float comparison (tests only) |

## Methodology

The reimplementation proceeded in two phases:

1. **Python teardown** -- module-by-module reverse engineering with stage-by-stage ground truth capture and validation against the Fortran binary. This produced the algorithmic understanding and the validation framework.

2. **Rust reimplementation** -- compiled implementation validated against the same ground truth. With algorithms understood from the Python phase, this became a translation problem rather than a discovery problem.

### Key Decisions

- **ndarray** for all 2D matrices (type-safe, column-major compatible)
- **Typed distance matrices** (Scale10/Scale100) preventing accidental precision mixing
- **f32 arithmetic** in PARSI for Fortran single-precision compatibility
- **Jacobi SVD** for 3x3 eigendecomposition (avoids LAPACK dependency)
- **Kabsch-Sander DSSP** from scratch (~420 lines) rather than wrapping external binary
- **Direct struct passing** between pipeline stages (no text serialization)

### Fortran Quirks Discovered

- **Stateful ca1 in DALICON**: `transrotate` modifies coordinates in-place; processing order matters
- **Fortran rounding**: `nint()` rounds half away from zero (not banker's rounding)
- **Single-precision intermediates**: float32 arithmetic despite float64 coordinates
- **FILTER95 stale state bug**: `lkeep` array retains values from previous protein when `idom > ndom`
- **PARSI score table sensitivity**: +/-1 int16 quantization differences propagate through branch-and-bound pruning

## License

This is a research project demonstrating AI-assisted software reimplementation methodology. The original DaliLite.v5 is by Liisa Holm (University of Helsinki).
