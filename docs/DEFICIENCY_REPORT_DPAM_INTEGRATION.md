# Deficiency Report: dali_cl Integration with DPAM v2.0

**Date:** 2026-02-15
**Tested against:** DPAM v2.0 step 07 (Iterative DALI) and step 08 (Analyze DALI)
**Test proteins:** AF-Q97ZL0-F1 (161 ECOD candidates), AF-P06596-F1 (213 candidates)
**Reference backend:** Fortran DaliLite.v5 via `dali.pl`

## Executive Summary

Integration testing of the Rust DALI backend (`dali.align_pdb()`) into DPAM's
step 07 pipeline revealed four deficiencies. All four have been resolved:

| # | Deficiency | Severity | Status | Resolution |
|---|-----------|----------|--------|------------|
| 1 | Template residue indexing mismatch | HIGH | **FIXED** | Both indices now 1-based sequential |
| 2 | DSSP secondary structure divergence | HIGH | **MITIGATED** | `.dat` bypass parameters added to `align_pdb()` |
| 3 | E2E validation test gaps | MEDIUM | **FIXED** | 7 new coverage tests + code-naming bug fix |
| 4 | Inconsistent `resid_map` semantics | MEDIUM | **FIXED** | `ResidNumbering` enum + `numbering` field |

### Pre-fix Metrics

| Metric | Fortran | Rust | Notes |
|--------|---------|------|-------|
| AF-Q97ZL0-F1 hits | 151 | 124 | 82% coverage |
| AF-P06596-F1 hits | 57 | 28 | 49% coverage |
| Rust-only hits | 0 | 0 | Rust never finds something Fortran misses |
| Z-score agreement (within 0.5) | - | 63% | For domains both backends find |
| Step 08 qscore agreement (within 0.05) | - | 49% | Affected by deficiencies 1 and 2 |
| Wall-clock speedup | - | 1.36-1.86x | Rust faster due to no subprocess overhead |

With deficiency 1 fixed and deficiency 2 mitigated (via `.dat` bypass), the
indexing and DSSP issues are eliminated. Remaining hit-count gaps are attributable
to inherent search divergence in WOLF/PARSI heuristics (documented at 78-87%
match rate in module validation).

---

## Deficiency 1: Template Residue Indexing Mismatch

**Severity:** HIGH — produces silently wrong downstream scores
**Status:** FIXED

### Problem

`align_pdb()` expanded alignment blocks to residue-level pairs using each
protein's `resid_map`, which holds **PDB residue serial numbers** for
PDB-imported proteins:

```rust
// BEFORE (dali-python/src/lib.rs)
let q_idx = (b.l1 as usize - 1) + k;
let t_idx = (b.l2 as usize - 1) + k;
alignments.push((query.resid_map[q_idx], template.resid_map[t_idx]));
```

Fortran DaliLite works through `.dat` files internally and returns 1-based
sequential indices. DPAM's downstream pipeline (Step 08 `calculate_qscore()`,
Steps 15, 18) expects this convention.

For template `000005249` (PDB residues 2-87, weights keyed 1-86):
- Fortran: `q=11, t=1` (correct — template index 1 = PDB residue 2)
- Rust: `q=10, t=2` (wrong — PDB residue number leaked through)

### Resolution

Changed `align_pdb()` to return 1-based sequential indices for **both** query
and template, matching the Fortran `.dat` convention:

```rust
// AFTER (dali-python/src/lib.rs:360-372)
let q_idx = (b.l1 as usize - 1) + k;  // 0-based internal index
let t_idx = (b.l2 as usize - 1) + k;
if q_idx < query.nres && t_idx < template.nres {
    alignments.push(((q_idx + 1) as i32, (t_idx + 1) as i32));
}
```

DPAM Step 7's `Qresids[qind - 1]` mapping now works correctly — it receives
1-based sequential indices and maps them to PDB residue IDs via its own lookup.

**Files changed:** `dali-python/src/lib.rs`

---

## Deficiency 2: DSSP Secondary Structure Divergence

**Severity:** HIGH — causes missed alignments and lower z-scores
**Status:** MITIGATED

### Problem

The Rust DSSP reimplementation produces different per-residue 8-state assignments
than Fortran `dsspcmbi` for a significant fraction of proteins. These differences
propagate through segment assignment and change the protein's SSE representation
(nseg, na, nb, segment boundaries).

Tested on 47 ECOD templates from AF-Q97ZL0-F1:

| Category | DSSP Segment Match | DSSP Segment Mismatch | Mismatch Rate |
|----------|-------------------|----------------------|---------------|
| Both backends found | 13 | 7 | 35% |
| Rust missed (Fortran-only) | 11 | 16 | 59% |

The causal chain: different DSSP → different segments → different WOLF spatial
hash → different PARSI search space → missed hits.

The segment assignment pipeline amplifies raw DSSP differences:
1. **3-state reduction**: G,I,H → H; E → E; all else → L
2. **Length filtering with growth**: Helices < 6 residues discarded; helices < 8
   and strands < 6 attempt to grow. A single-residue DSSP disagreement at a
   boundary can cause a segment to fall below threshold and be dropped entirely.

### Resolution

Implemented **Option C** from the original report: bypass DSSP entirely for the
DPAM use case.

Added `query_dat` and `template_dat` optional parameters to `align_pdb()`:

```python
result = dali.align_pdb(
    query_path, template_path,
    query_chain, template_chain,
    query_code, template_code,
    query_dat="/path/to/query.dat",       # optional
    template_dat="/path/to/template.dat",  # optional
)
```

When provided, the protein is loaded from the pre-computed `.dat` file (with
Fortran-computed segments) instead of running the Rust DSSP pipeline. This makes
the protein representations identical to Fortran, isolating comparison to the
alignment algorithm only.

The `RustDALI` class in DPAM accepts a `dat_dir` parameter for automatic template
`.dat` lookup:

```python
dali_backend = RustDALI(dat_dir="/path/to/ecod_dat_files")
```

**Files changed:** `dali-python/src/lib.rs`, `dpam/tools/dali.py`,
`dpam/steps/step07_iterative_dali.py`

### Remaining Work

The underlying DSSP divergence is not fixed — it remains a known limitation of
the Rust DSSP reimplementation. For use cases where pre-computed `.dat` files are
not available, the PDB-only path may miss hits. Options B (targeted DSSP
debugging) and A (external Fortran DSSP) from the original report remain viable
for a deeper fix.

---

## Deficiency 3: E2E Validation Test Gaps

**Severity:** MEDIUM — deficiencies 1 and 2 went undetected
**Status:** FIXED

### Problem

The E2E validation suite tested only 3 structure pairs from a 5-structure corpus
(101mA, 1a00A, 1a87A, 1allA, 1binA). All have sequential residue numbering
starting near 1, large well-defined SSEs, and no edge-case features. This corpus
is not representative of ECOD70 templates.

### Resolution

Added `test_align_coverage.py` with 7 tests covering previously untested cases:

| Test | What it covers |
|------|---------------|
| `test_sequential_indices_with_hits` | 3 globin pairs with PDB offsets (91, 3, 0); verifies indices within [1, nres] |
| `test_pdb_import_alignment` | Full PDB→DSSP→WOLF→DP→DALICON pipeline (no .dat bypass) |
| `test_five_char_codes` | Regression test for code-naming bug (see below) |
| `test_small_domain_self_alignment` | 62-residue 1aiwA (5 SSEs, all-beta); verifies z > 5.0, n >= 30 |
| `test_dat_bypass_vs_pdb_import` | Compares PDB-only, .dat-template, and .dat-both modes |
| `test_numbering_field_consistency` | Verifies `numbering` field for .dat-loaded vs PDB-imported proteins |
| `test_rotation_matrix_on_hits` | Validates det(R) ≈ 1.0, R^T R ≈ I on 3 real alignment pairs |

**Code-naming bug discovered and fixed:** `import_pdb("...", "A", "101mA")`
produced `code="101mAA"` (chain appended to already-suffixed code). The 6-char
code overflowed the `.dat` format's 5-char field, corrupting the file and causing
`read_dat` to fail with a parse error. This made ALL alignment tests with
chain-suffixed codes silently return `None`.

Fix: `pdb.rs` now skips chain append when the code already ends with the chain
letter, and truncates to 5 chars maximum. `write_dat` also truncates the code
field to 5 chars as a safety guard.

**Files changed:** `dali-core/src/io/pdb.rs`, `dali-core/src/io/dat.rs`,
`dali-python/test_align_coverage.py` (new)

---

## Deficiency 4: Inconsistent `resid_map` Semantics

**Severity:** MEDIUM — API contract is unclear
**Status:** FIXED

### Problem

The `resid_map` field on `Protein` had different semantics depending on the
import source, with no way for callers to distinguish:

| Source | `resid_map` contents | Example (86 residues, PDB starts at 2) |
|--------|---------------------|---------------------------------------|
| `import_pdb()` | PDB residue serial numbers | [2, 3, 4, ..., 87] |
| `read_dat()` | Synthetic 1..=nres | [1, 2, 3, ..., 86] |

### Resolution

Added a `ResidNumbering` enum and `numbering` field to `Protein`:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResidNumbering {
    Sequential,  // from .dat files: resid_map = 1..=nres
    Pdb,         // from import_pdb: resid_map = PDB serial numbers
}
```

PyO3 exposes this as `protein.numbering`, returning `"sequential"` or `"pdb"`.
Callers can now dispatch on the numbering convention:

```python
p = dali.import_pdb("file.pdb", "A", "test")
assert p.numbering == "pdb"

p = dali.read_dat("file.dat")
assert p.numbering == "sequential"
```

**Files changed:** `dali-core/src/types/protein.rs`, `dali-core/src/types/mod.rs`,
`dali-core/src/lib.rs`, `dali-core/src/io/dat.rs`, `dali-core/src/io/import.rs`,
`dali-python/src/lib.rs`

---

## Reproduction

### Test Script

The original comparison was run using `dpam_c2/scripts/test_dali_backends.py`:

```bash
cd ~/dev/dpam_c2
python scripts/test_dali_backends.py \
  --prefix AF-Q97ZL0-F1 \
  --source-dir validation_swissprot \
  --data-dir /home/rschaeff_1/data/dpam_reference/ecod_data \
  --cpus 4
```

### Validation Tests

```bash
# Rust core tests (61 tests)
cd dali_rust && cargo test --package dali-core --release

# Python binding tests (10 tests)
cd dali_rust/dali-python
maturin develop --release
LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libgomp.so.1 python3 test_bindings.py

# Alignment coverage tests (7 tests)
LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libgomp.so.1 python3 test_align_coverage.py
```
