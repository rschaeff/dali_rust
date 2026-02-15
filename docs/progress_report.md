# DaliLite.v5 Reimplementation — Progress Report

**Date:** 2026-02-13
**Project:** AI-assisted validated reimplementation of DaliLite.v5 structural comparison pipeline

## Executive Summary

Five of six pipeline modules have been reimplemented in Python and validated against
captured Fortran ground truth. The reimplementation totals ~6,200 lines of Python
(from ~8,100 lines of Fortran 77), plus ~2,200 lines of test and validation framework.

| Module    | Python Lines | Fortran Lines | Validation (5-struct) | Validation (18-struct) |
|-----------|-------------|---------------|----------------------|----------------------|
| WOLF      | 1,186       | ~860          | 21/21 (100%)         | 279/279 (100%)       |
| DP        | 559         | ~600          | 10/10 (100%)         | 24/24 (100%)         |
| DALICON   | 1,507       | ~800          | 10/10 (100%)         | 24/24 (100%)         |
| PARSI     | 2,400       | ~4,765        | 379/435 (87%)        | not yet              |
| FILTER95  | 538         | ~600          | 38/38 (100%)         | not yet              |
| **Total** | **6,190**   | **~7,625**    |                      |                      |

Validation + test framework: 2,177 lines.
Grand total Python: **8,367 lines**.

## Pipeline Coverage

The DaliLite structural comparison pipeline:

```
Input (.dat) → WOLF → DP → DALICON → PARSI → FILTER95 → final DALICON → DP → Output
```

**Completed modules:**
1. **WOLF** — SSE hashing and spatial voting for initial alignment seed (7 files)
2. **DP** — Z-score computation via dynamic programming (2 files)
3. **DALICON** — Genetic algorithm alignment refinement (5 files)
4. **PARSI** — Exhaustive branch-and-bound alignment optimization (8 files)
5. **FILTER95** — Redundancy filtering with FITZ superposition refinement (2 files)

**Remaining:**
6. Final DALICON → DP pass (reuses existing modules, just different parameters)
7. SOAP (NW on CA distances for proteins with <3 SSEs — rare edge case)
8. DSSP/import pipeline (not priority — using existing .dat files)

## Validation Framework

Ground truth is captured stage-by-stage from the reference DaliLite.v5 Fortran binary
using `validation/capture.py`. Each module's output is compared against the captured
reference using tolerance-based matching:

- **WOLF**: Exact integer match on nblock and ranges
- **DP**: Score tolerance of 2.0, Z-score tolerance
- **DALICON**: Score and Z-score tolerance matching
- **PARSI**: Score tolerance of 100 (suboptimal alignments differ due to stack ordering)
- **FILTER95**: Score tolerance of 1% or 100 units, Z-score tolerance of 0.5

### Test Corpus

**Original (5 structures):** 101mA, 1a00A, 1a87A, 1allA, 1binA
- Size range: 141–297 residues
- 25 directed pairs per module
- Ground truth for all 6 pipeline stages

**Expanded (18 structures, ECOD-diverse):**
- Selected from ECOD database spanning all major fold types (a.1 through a.18)
- Size range: 62–456 residues, 4–30 SSEs
- Ground truth captured for WOLF → DP → DALICON stages

## Key Implementation Challenges

### WOLF (solved)
- Fortran 1-based indexing for .dat file parsing
- `nint()` rounding differences (Fortran rounds 0.5 up; Python uses banker's rounding)
- Grid coordinate clamping for spatial hash

### DP (solved)
- Score tolerance needed for floating-point arithmetic differences

### DALICON (solved)
- Complex genetic algorithm with multiple scoring stages
- `testi` evaluation function with cross-validation-like fold scoring

### PARSI (87% validated, root cause identified)
- **4 critical 0-based vs 1-based bugs** found and fixed during development:
  1. `s_beg`/`s_end` offset (+29) missing in singletex/segsegscore
  2. `dostack` lprint parameter — duplicate output from 10-block and perresiduescore paths
  3. `compress` left-shift off-by-one (1-based vs 0-based array offsets)
  4. `split` sentinel value — segment 0 is valid in Python (changed sentinel from 0 to -1)
- 56/435 refine lines unmatched at tolerance=100 (431/435 = 99.1% at tolerance=500K)
- **Root cause: heuristic branch-and-bound search divergence**, NOT implementation bugs
  - Score lookup table has 194/64,000 entries that differ by ±1 (float32 arithmetic)
  - These ±1 differences change priority stack ordering in the branch-and-bound
  - Different search paths → different suboptimal alignments found
  - Python sometimes finds genuinely higher-scoring alignments (independently verified)
  - 93% of mismatches involve 1a87A (297 residues, most complex domain tree)
- See `docs/parsi_mismatch_analysis.md` for full analysis
- Performance: ~65-80 min for 25 pairs in Python vs seconds in Fortran

### FILTER95 (solved)
- **Critical discovery:** In the Fortran module version (comparemodules.f), `oldidom` is a
  local variable initialized to 0 on every subroutine call. Since the subroutine is called
  once per input line, this means:
  - The dedup check never triggers (idom is never 0 for valid domains)
  - Domain-local data (resix/xiser/xca) is always rebuilt from scratch
  - Multiple output lines per (cd1, cd2, idom) triplet are expected
- `lnew1` flag must persist across skipped (non-kept domain) lines
- Custom `fitz95` function needed to match Fortran's initial-alignment-aware convergence check
- Fast-path Z-score bypass is COMMENTED OUT in the module version — all alignments go through FITZ

## Architecture

```
dali_cl/
├── wolf/           # 7 files, 1186 lines — SSE hashing + spatial voting
│   ├── __init__.py
│   ├── dat_reader.py      # .dat file parser → Protein dataclass
│   ├── geometry.py        # SSE midpoints, directions, canonical frames
│   ├── spatial_hash.py    # 3D grid for rapid descriptor lookup
│   ├── alignment.py       # NW dynamic programming, fitz iterative refinement
│   ├── superposition.py   # Kabsch (u3b), transrotate
│   └── wolf.py            # Top-level orchestration
├── dp/             # 2 files, 559 lines — Z-score computation
│   ├── __init__.py
│   └── dp.py
├── dalicon/        # 5 files, 1507 lines — Genetic algorithm refinement
│   ├── __init__.py
│   ├── dalicon.py         # Entry point, orchestration
│   ├── lean.py            # GA core: crossover, mutation, selection
│   ├── scoring.py         # Scoring functions
│   ├── output.py          # DCCP output formatting
│   └── testi.py           # Evaluation function
├── parsi/          # 8 files, 2400 lines — Branch-and-bound optimization
│   ├── __init__.py
│   ├── parsi.py           # Entry points: dowork_parsi, run_parsi
│   ├── align.py           # Core: align(), dostack(), compress(), _write_refine
│   ├── scoring.py         # Score tables: singletex, segsegscore, doubletex
│   ├── stack.py           # PriorityStack, getnextbest, declump
│   ├── domain_tree.py     # Domain decomposition
│   ├── protein_io.py      # Protein data I/O, distance matrices
│   └── parsizes.py        # Constants
├── filter95/       # 2 files, 538 lines — Redundancy filtering + FITZ
│   ├── __init__.py
│   └── filter95.py        # All functions: readproteindata95, scorefun95, fitz95, etc.
├── validation/     # 3 files, 1252 lines — Ground truth framework
│   ├── capture.py         # Capture ground truth from DaliLite binary
│   ├── compare.py         # Compare Python output to ground truth
│   ├── parsers.py         # Output format parsers
│   └── fixtures/          # Captured ground truth data
├── test_wolf.py    # WOLF validation test
├── test_dp.py      # DP validation test
├── test_dalicon.py # DALICON validation test
├── test_parsi.py   # PARSI validation test
├── test_filter95.py        # FILTER95 validation test
└── test_expanded_corpus.py # 18-structure expanded corpus test
```

## Next Steps

1. **Investigate PARSI 56 missing entries** — concentrated in 1a87A (297-residue) pairs;
   likely float32 vs float64 divergence in inner scoring loops
2. **Scale PARSI to 18-structure corpus** — requires capturing ground truth
3. **Scale FILTER95 to 18-structure corpus** — requires PARSI ground truth first
4. **Final DALICON → DP pass** — reuses existing modules with different parameters
5. **End-to-end pipeline integration** — chain all modules for a complete comparison run
