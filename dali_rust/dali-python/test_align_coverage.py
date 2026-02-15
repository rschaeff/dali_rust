"""
Expanded alignment coverage tests (Deficiency 3 fix).

Tests align_pdb() on structures that produce actual alignments, covering
cases the original 5-structure E2E tests missed:

- PDB numbering not starting at 1 (1a87A starts at 91, 1allA at 3, 101mA at 0)
- 1-based sequential indices verified against nres bounds
- .dat bypass vs PDB import (DSSP divergence detection)
- Numbering field consistency across import paths
- Rotation matrix validity on real alignments
- Small domain alignment (1aiwA: 62 res, 5 SSE self-comparison)
- 5-char code handling (import_pdb with chain-suffixed codes)

NOTE: Some globin pairs (101mA/1binA, 1a87A/1allA) only produce hits when using
.dat files because Rust DSSP assigns different segments than Fortran dsspcmbi.
Tests use .dat bypass where needed to ensure alignments are exercised.
"""

import os
import sys
import traceback
import numpy as np

import dali

# Paths
ORIG_PDB = os.path.join(os.path.dirname(__file__),
    "../../validation/fixtures/structures/pdb")
ORIG_DAT = os.path.join(os.path.dirname(__file__),
    "../../validation/fixtures/structures")
CORPUS_PDB = os.path.join(os.path.dirname(__file__),
    "../../validation/corpus_expansion/pdb")
CORPUS_DAT = os.path.join(os.path.dirname(__file__),
    "../../validation/corpus_expansion/ground_truth_18/structures")

# Original 5-structure corpus (globins — produce known alignments)
# All confirmed E2E pairs: 101mA/1a00A, 101mA/1binA, 1a87A/1allA
ORIG_5 = [
    ("pdb101m.ent.gz", "A", "101mA"),   # nres=154, starts at 0
    ("pdb1a00.ent.gz", "A", "1a00A"),    # nres=141, starts at 1
    ("pdb1a87.ent.gz", "A", "1a87A"),    # nres=297, starts at 91 (big offset)
    ("pdb1all.ent.gz", "A", "1allA"),    # nres=160, starts at 3
    ("pdb1bin.ent.gz", "A", "1binA"),    # nres=143, starts at 1
]


def _dat_path(code):
    return os.path.join(ORIG_DAT, f"{code}.dat")


def test_sequential_indices_with_hits():
    """align_pdb returns 1-based sequential indices for pairs that produce hits.

    Uses globin pairs with .dat bypass (to ensure hits regardless of DSSP
    divergence), and verifies indices are within [1, nres] — not PDB numbers.
    """
    test_pairs = [
        ("1a87A", "1allA"),   # query starts at 91, template at 3
        ("101mA", "1a00A"),   # query starts at 0, template at 1
        ("101mA", "1binA"),   # query starts at 0, template at 1
    ]

    hits_found = 0
    for q_code, t_code in test_pairs:
        q_entry = next(e for e in ORIG_5 if e[2] == q_code)
        t_entry = next(e for e in ORIG_5 if e[2] == t_code)
        q_path = os.path.join(ORIG_PDB, q_entry[0])
        t_path = os.path.join(ORIG_PDB, t_entry[0])

        if not (os.path.exists(q_path) and os.path.exists(t_path)):
            print(f"    SKIPPED {q_code} vs {t_code}: files not found")
            continue

        # Use .dat bypass for both to ensure hits (DSSP divergence blocks some PDB-only pairs)
        result = dali.align_pdb(
            q_path, t_path, q_entry[1], t_entry[1], q_code, t_code,
            query_dat=_dat_path(q_code), template_dat=_dat_path(t_code))

        assert result is not None, \
            f"{q_code} vs {t_code}: expected alignment (known globin pair with .dat), got None"

        # Load proteins to check nres (from .dat for consistent nres)
        q_prot = dali.read_dat(_dat_path(q_code))
        t_prot = dali.read_dat(_dat_path(t_code))

        # Verify all indices are 1-based sequential (within nres)
        for q_idx, t_idx in result.alignments:
            assert 1 <= q_idx <= q_prot.nres, \
                f"{q_code} vs {t_code}: query index {q_idx} out of range [1, {q_prot.nres}]"
            assert 1 <= t_idx <= t_prot.nres, \
                f"{q_code} vs {t_code}: template index {t_idx} out of range [1, {t_prot.nres}]"

        q_max_idx = max(q for q, t in result.alignments)
        t_max_idx = max(t for q, t in result.alignments)
        assert q_max_idx <= q_prot.nres, \
            f"{q_code}: max query index {q_max_idx} > nres {q_prot.nres} (PDB numbering leaked)"
        assert t_max_idx <= t_prot.nres, \
            f"{t_code}: max template index {t_max_idx} > nres {t_prot.nres} (PDB numbering leaked)"

        print(f"    {q_code} vs {t_code}: z={result.zscore:.1f}, n={result.n_aligned}, "
              f"q_range=[1..{q_max_idx}]/{q_prot.nres}, t_range=[1..{t_max_idx}]/{t_prot.nres}")
        hits_found += 1

    assert hits_found == len(test_pairs), f"only {hits_found}/{len(test_pairs)} pairs produced hits"


def test_pdb_import_alignment():
    """align_pdb works PDB-only (no .dat bypass) for pairs where Rust DSSP succeeds.

    101mA vs 1a00A reliably produces hits with PDB import (nseg=7/6 vs dat 5/5).
    This validates the full PDB→DSSP→segments→WOLF→DP→DALICON pipeline.
    """
    q_path = os.path.join(ORIG_PDB, "pdb101m.ent.gz")
    t_path = os.path.join(ORIG_PDB, "pdb1a00.ent.gz")
    if not (os.path.exists(q_path) and os.path.exists(t_path)):
        print("    SKIPPED: PDB files not found")
        return

    result = dali.align_pdb(q_path, t_path, "A", "A", "101mA", "1a00A")
    assert result is not None, "101mA vs 1a00A (PDB import) returned None"
    assert result.zscore > 10.0, f"expected z > 10.0 for globin pair, got {result.zscore:.1f}"
    assert result.n_aligned > 100, f"expected n > 100, got {result.n_aligned}"

    # Indices should be 1-based sequential
    q_prot = dali.import_pdb(q_path, "A", "101mA")
    for q_idx, t_idx in result.alignments:
        assert 1 <= q_idx <= q_prot.nres, f"query index {q_idx} out of range"

    print(f"    101mA vs 1a00A (PDB): z={result.zscore:.1f}, n={result.n_aligned}")


def test_five_char_codes():
    """align_pdb works correctly with 5-char codes (chain-suffixed).

    Regression test for bug where import_pdb("...", "A", "101mA") produced
    code="101mAA" (chain appended twice), overflowing the .dat 5-char field.
    """
    q_path = os.path.join(ORIG_PDB, "pdb101m.ent.gz")
    t_path = os.path.join(ORIG_PDB, "pdb1a00.ent.gz")
    if not (os.path.exists(q_path) and os.path.exists(t_path)):
        print("    SKIPPED: PDB files not found")
        return

    # 5-char codes (chain-suffixed)
    result_5 = dali.align_pdb(q_path, t_path, "A", "A", "101mA", "1a00A")
    assert result_5 is not None, "5-char codes returned None"

    # 4-char codes (chain appended by import_pdb)
    result_4 = dali.align_pdb(q_path, t_path, "A", "A", "101m", "1a00")
    assert result_4 is not None, "4-char codes returned None"

    # Both should produce identical results
    assert abs(result_5.zscore - result_4.zscore) < 0.1, \
        f"z-score mismatch: 5-char={result_5.zscore:.1f} vs 4-char={result_4.zscore:.1f}"
    assert result_5.n_aligned == result_4.n_aligned, \
        f"n_aligned mismatch: 5-char={result_5.n_aligned} vs 4-char={result_4.n_aligned}"

    # Verify import_pdb code construction
    p5 = dali.import_pdb(q_path, "A", "101mA")
    p4 = dali.import_pdb(q_path, "A", "101m")
    assert p5.code == "101mA", f"5-char import produced code={p5.code!r}"
    assert p4.code == "101mA", f"4-char import produced code={p4.code!r}"

    print(f"    5-char: z={result_5.zscore:.1f}, n={result_5.n_aligned}")
    print(f"    4-char: z={result_4.zscore:.1f}, n={result_4.n_aligned}")
    print(f"    codes: import('101mA','A')={p5.code!r}, import('101m','A')={p4.code!r}")


def test_small_domain_self_alignment():
    """align_pdb works for small domains via self-comparison with different codes."""
    aiw_path = os.path.join(CORPUS_PDB, "pdb1aiw.ent.gz")
    if not os.path.exists(aiw_path):
        print("    SKIPPED: 1aiwA PDB not found")
        return

    p = dali.import_pdb(aiw_path, "A", "1aiwA")
    print(f"    1aiwA: nres={p.nres}, nseg={p.nseg}, na={p.na}, nb={p.nb}")
    assert p.nseg >= 3, f"nseg={p.nseg} too low for WOLF (need >= 3)"

    # Self-comparison needs different codes (compare_pair filters cd1==cd2).
    # Use 4-char codes so import_pdb appends chain → 5-char codes "querA" / "tmplA".
    result = dali.align_pdb(aiw_path, aiw_path, "A", "A", "quer", "tmpl")
    assert result is not None, "1aiwA self-alignment returned None"
    assert result.n_aligned >= 30, \
        f"self-alignment of 62-residue protein should align >= 30, got {result.n_aligned}"
    assert result.zscore > 5.0, \
        f"self-alignment should have high z-score, got {result.zscore:.1f}"

    # Verify indices are within bounds
    for q_idx, t_idx in result.alignments:
        assert 1 <= q_idx <= p.nres and 1 <= t_idx <= p.nres, \
            f"index out of range: ({q_idx}, {t_idx}) for nres={p.nres}"

    print(f"    1aiwA self-align: z={result.zscore:.1f}, n={result.n_aligned}/{p.nres}")


def test_dat_bypass_vs_pdb_import():
    """Compare align_pdb results with and without template .dat bypass.

    Uses 101mA vs 1a00A (the one pair that works PDB-only). Both paths should
    produce valid alignments; z-scores may differ if DSSP differs.
    """
    q_path = os.path.join(ORIG_PDB, "pdb101m.ent.gz")
    t_pdb = os.path.join(ORIG_PDB, "pdb1a00.ent.gz")
    t_dat = os.path.join(ORIG_DAT, "1a00A.dat")

    if not (os.path.exists(q_path) and os.path.exists(t_pdb) and os.path.exists(t_dat)):
        print("    SKIPPED: files not found")
        return

    # PDB-only (Rust DSSP for both)
    result_pdb = dali.align_pdb(q_path, t_pdb, "A", "A", "101mA", "1a00A")

    # .dat bypass for template only
    result_dat = dali.align_pdb(
        q_path, t_pdb, "A", "A", "101mA", "1a00A", template_dat=t_dat)

    # .dat bypass for both
    q_dat = os.path.join(ORIG_DAT, "101mA.dat")
    result_both = dali.align_pdb(
        q_path, t_pdb, "A", "A", "101mA", "1a00A",
        query_dat=q_dat, template_dat=t_dat)

    assert result_pdb is not None, "PDB-only alignment returned None"
    assert result_dat is not None, ".dat-template alignment returned None"
    assert result_both is not None, ".dat-both alignment returned None"

    # All should produce valid alignments with positive z-scores
    for label, r in [("pdb", result_pdb), ("dat-template", result_dat), ("dat-both", result_both)]:
        assert r.zscore > 0, f"{label}: zscore <= 0"
        for q, t in r.alignments:
            assert q >= 1 and t >= 1, f"{label}: index < 1"

    print(f"    pdb:          z={result_pdb.zscore:.1f} n={result_pdb.n_aligned}")
    print(f"    dat-template: z={result_dat.zscore:.1f} n={result_dat.n_aligned}")
    print(f"    dat-both:     z={result_both.zscore:.1f} n={result_both.n_aligned}")


def test_numbering_field_consistency():
    """Verify numbering field is set correctly for all import paths."""
    # .dat-loaded: should be "sequential"
    for code in ["101mA", "1a00A", "1binA"]:
        dat_path = os.path.join(ORIG_DAT, f"{code}.dat")
        if not os.path.exists(dat_path):
            continue
        p = dali.read_dat(dat_path)
        assert p.numbering == "sequential", \
            f"{code} from .dat: expected 'sequential', got '{p.numbering}'"
        assert p.resid_map == list(range(1, p.nres + 1)), \
            f"{code} from .dat: resid_map not sequential"

    # PDB-imported: should be "pdb"
    for pdb_file, chain, code in [
        ("pdb1a87.ent.gz", "A", "1a87A"),  # starts at 91
        ("pdb1all.ent.gz", "A", "1allA"),   # starts at 3
        ("pdb101m.ent.gz", "A", "101mA"),   # starts at 0
    ]:
        pdb_path = os.path.join(ORIG_PDB, pdb_file)
        if not os.path.exists(pdb_path):
            continue
        p = dali.import_pdb(pdb_path, chain, code)
        assert p.numbering == "pdb", \
            f"{code} from PDB: expected 'pdb', got '{p.numbering}'"
        print(f"    {code}: numbering={p.numbering}, resid_map[0]={p.resid_map[0]}, "
              f"resid_map[-1]={p.resid_map[-1]}")

    # .dat from 18-corpus
    for code in ["1a7sA", "1bbhA", "1bcoA"]:
        dat_path = os.path.join(CORPUS_DAT, f"{code}.dat")
        if not os.path.exists(dat_path):
            continue
        p = dali.read_dat(dat_path)
        assert p.numbering == "sequential", \
            f"{code} corpus .dat: expected 'sequential', got '{p.numbering}'"


def test_rotation_matrix_on_hits():
    """Verify rotation matrices are valid proper rotations for real alignments."""
    pairs = [
        ("101mA", "1a00A"),   # globin pair
        ("101mA", "1binA"),   # globin pair
        ("1a87A", "1allA"),   # globin pair (with PDB offset)
    ]

    validated = 0
    for q_code, t_code in pairs:
        q_entry = next(e for e in ORIG_5 if e[2] == q_code)
        t_entry = next(e for e in ORIG_5 if e[2] == t_code)
        q_path = os.path.join(ORIG_PDB, q_entry[0])
        t_path = os.path.join(ORIG_PDB, t_entry[0])

        if not (os.path.exists(q_path) and os.path.exists(t_path)):
            continue

        # Use .dat bypass to ensure hits
        result = dali.align_pdb(
            q_path, t_path, q_entry[1], t_entry[1], q_code, t_code,
            query_dat=_dat_path(q_code), template_dat=_dat_path(t_code))
        assert result is not None, f"{q_code} vs {t_code}: expected hit"

        R = np.array(result.rotation)
        det = np.linalg.det(R)
        assert abs(det - 1.0) < 0.01, \
            f"{q_code} vs {t_code}: rotation det={det:.4f}, expected ~1.0"

        eye = R.T @ R
        ident_err = np.abs(eye - np.eye(3)).max()
        assert ident_err < 0.01, \
            f"{q_code} vs {t_code}: R^T R not identity, max err={ident_err:.4f}"

        assert len(result.translation) == 3
        assert all(np.isfinite(result.translation))

        print(f"    {q_code} vs {t_code}: z={result.zscore:.1f}, "
              f"n={result.n_aligned}, det={det:.6f}")
        validated += 1

    assert validated == len(pairs), f"only {validated}/{len(pairs)} pairs validated"


if __name__ == "__main__":
    tests = [
        ("sequential_indices_with_hits",
         test_sequential_indices_with_hits),
        ("pdb_import_alignment",
         test_pdb_import_alignment),
        ("five_char_codes",
         test_five_char_codes),
        ("small_domain_self_alignment",
         test_small_domain_self_alignment),
        ("dat_bypass_vs_pdb_import",
         test_dat_bypass_vs_pdb_import),
        ("numbering_field_consistency",
         test_numbering_field_consistency),
        ("rotation_matrix_on_hits",
         test_rotation_matrix_on_hits),
    ]

    passed = 0
    failed = 0
    for name, fn in tests:
        try:
            fn()
            print(f"PASS: {name}")
            passed += 1
        except Exception as e:
            traceback.print_exc()
            print(f"FAIL: {name}")
            failed += 1

    print(f"\n{passed}/{passed + failed} tests passed")
    if failed:
        sys.exit(1)
