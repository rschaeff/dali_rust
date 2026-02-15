"""Validation tests for dali PyO3 bindings."""

import os
import sys
import traceback
import numpy as np

import dali

DAT_DIR = os.path.join(os.path.dirname(__file__), "../../validation/fixtures/structures")
PDB_PATH = os.path.join(DAT_DIR, "pdb/pdb101m.ent.gz")


def test_read_dat():
    """read_dat returns Protein with correct fields."""
    p = dali.read_dat(os.path.join(DAT_DIR, "101mA.dat"))
    assert p.code == "101mA", f"expected code '101mA', got '{p.code}'"
    assert p.nres == 154, f"expected nres 154, got {p.nres}"
    assert p.nseg > 0, f"expected nseg > 0, got {p.nseg}"
    print(f"  read_dat: {p}")


def test_ca_coords():
    """ca_coords() returns numpy array with shape (3, nres)."""
    p = dali.read_dat(os.path.join(DAT_DIR, "101mA.dat"))
    ca = p.ca_coords()
    assert isinstance(ca, np.ndarray), f"expected ndarray, got {type(ca)}"
    assert ca.shape == (3, p.nres), f"expected shape (3, {p.nres}), got {ca.shape}"
    assert ca.dtype == np.float64
    assert np.all(np.isfinite(ca))
    print(f"  ca_coords: shape={ca.shape}, dtype={ca.dtype}")


def test_import_pdb():
    """import_pdb returns Protein from PDB file."""
    if not os.path.exists(PDB_PATH):
        print(f"  import_pdb: SKIPPED (no PDB file at {PDB_PATH})")
        return
    p = dali.import_pdb(PDB_PATH, "A", "101m")
    assert p.nres > 0, f"expected nres > 0, got {p.nres}"
    assert len(p.sequence) == p.nres, f"expected seq len {p.nres}, got {len(p.sequence)}"
    print(f"  import_pdb: {p}")


def test_protein_store():
    """ProteinStore loads .dat files on demand."""
    store = dali.ProteinStore(DAT_DIR)
    assert "101mA" in store, "expected '101mA' in store"
    assert "XXXXX" not in store, "expected 'XXXXX' not in store"
    p = store.get_protein("101mA")
    assert p.nres == 154, f"expected nres 154, got {p.nres}"
    assert len(store) >= 1, f"expected len >= 1, got {len(store)}"
    print(f"  ProteinStore: {store}")


def test_compare_pair():
    """compare_pair produces results matching E2E reference."""
    store = dali.ProteinStore(DAT_DIR)
    results = dali.compare_pair("101mA", "1a00A", store)
    assert len(results) > 0, "expected at least one result"
    assert len(results) == 4, f"expected 4 results, got {len(results)}"
    for r in results:
        assert r.score > 0
        assert r.nblock > 0
        assert len(r.blocks) == r.nblock
    print(f"  compare_pair: {len(results)} results")
    for r in results:
        print(f"    {r}")


def test_run_wolf_path():
    """run_wolf_path produces results."""
    store = dali.ProteinStore(DAT_DIR)
    results = dali.run_wolf_path("101mA", ["1a00A", "1binA"], store)
    assert len(results) > 0, "expected at least one wolf path result"
    print(f"  run_wolf_path: {len(results)} results")
    for r in results:
        print(f"    {r}")


def test_resid_map():
    """resid_map returns sequential IDs for .dat-loaded, PDB IDs for imported."""
    # .dat-loaded: sequential 1..=nres, numbering="sequential"
    p_dat = dali.read_dat(os.path.join(DAT_DIR, "101mA.dat"))
    rmap = p_dat.resid_map
    assert len(rmap) == p_dat.nres, f"expected len {p_dat.nres}, got {len(rmap)}"
    assert rmap == list(range(1, p_dat.nres + 1)), "expected sequential for .dat-loaded"
    assert p_dat.numbering == "sequential", f"expected 'sequential', got '{p_dat.numbering}'"

    # PDB-imported: should be PDB residue numbers, numbering="pdb"
    if os.path.exists(PDB_PATH):
        p_pdb = dali.import_pdb(PDB_PATH, "A", "101m")
        rmap_pdb = p_pdb.resid_map
        assert len(rmap_pdb) == p_pdb.nres
        assert all(isinstance(x, int) for x in rmap_pdb)
        assert p_pdb.numbering == "pdb", f"expected 'pdb', got '{p_pdb.numbering}'"
        print(f"  resid_map: dat={rmap[:3]}... ({p_dat.numbering})  pdb={rmap_pdb[:3]}... ({p_pdb.numbering})")
    else:
        print(f"  resid_map: dat={rmap[:3]}... ({p_dat.numbering}) (PDB skipped)")


def test_write_dat_roundtrip():
    """write_dat produces a file that read_dat can parse back correctly."""
    import tempfile
    p_orig = dali.read_dat(os.path.join(DAT_DIR, "101mA.dat"))
    with tempfile.NamedTemporaryFile(suffix=".dat", delete=False) as f:
        tmp_path = f.name
    try:
        p_orig.write_dat(tmp_path)
        p_back = dali.read_dat(tmp_path)
        assert p_back.code == p_orig.code, f"code mismatch: {p_back.code} vs {p_orig.code}"
        assert p_back.nres == p_orig.nres, f"nres mismatch: {p_back.nres} vs {p_orig.nres}"
        assert p_back.nseg == p_orig.nseg, f"nseg mismatch: {p_back.nseg} vs {p_orig.nseg}"
        assert p_back.na == p_orig.na, f"na mismatch: {p_back.na} vs {p_orig.na}"
        assert p_back.nb == p_orig.nb, f"nb mismatch: {p_back.nb} vs {p_orig.nb}"
        # CA coords
        ca_orig = p_orig.ca_coords()
        ca_back = p_back.ca_coords()
        max_diff = np.abs(ca_orig - ca_back).max()
        assert max_diff < 0.1, f"CA coord max diff {max_diff:.3f} >= 0.1"
        print(f"  write_dat roundtrip: OK (max CA diff={max_diff:.4f})")
    finally:
        os.unlink(tmp_path)


def test_align_pdb():
    """align_pdb returns AlignResult with rotation, translation, alignments."""
    pdb_dir = os.path.join(DAT_DIR, "pdb")
    q_path = os.path.join(pdb_dir, "pdb101m.ent.gz")
    t_path = os.path.join(pdb_dir, "pdb1a00.ent.gz")
    if not (os.path.exists(q_path) and os.path.exists(t_path)):
        print(f"  align_pdb: SKIPPED (PDB files not found)")
        return

    result = dali.align_pdb(q_path, t_path, "A", "A", "101m", "1a00")
    assert result is not None, "expected AlignResult, got None"
    assert result.zscore > 0, f"expected zscore > 0, got {result.zscore}"
    assert result.n_aligned > 0, f"expected n_aligned > 0, got {result.n_aligned}"
    assert len(result.alignments) == result.n_aligned
    assert len(result.rotation) == 3
    assert all(len(row) == 3 for row in result.rotation)
    assert len(result.translation) == 3
    assert len(result.blocks) > 0

    # Verify rotation is a proper rotation matrix (det ~= 1, R^T R ~= I)
    R = np.array(result.rotation)
    det = np.linalg.det(R)
    assert abs(det - 1.0) < 0.01, f"rotation det={det:.4f}, expected ~1.0"
    eye = R.T @ R
    ident_err = np.abs(eye - np.eye(3)).max()
    assert ident_err < 0.01, f"R^T R not identity, max err={ident_err:.4f}"

    # Alignments should be 1-based sequential indices (matching Fortran .dat convention)
    for q_res, t_res in result.alignments:
        assert isinstance(q_res, int) and isinstance(t_res, int)
        assert q_res >= 1, f"query index {q_res} < 1"
        assert t_res >= 1, f"template index {t_res} < 1"

    print(f"  align_pdb: {result}")
    print(f"    n_aligned={result.n_aligned}, nblock={len(result.blocks)}")
    print(f"    first 3 pairs: {result.alignments[:3]}")


def test_align_pdb_with_dat():
    """align_pdb with template_dat bypasses DSSP and uses .dat file directly."""
    pdb_dir = os.path.join(DAT_DIR, "pdb")
    q_path = os.path.join(pdb_dir, "pdb101m.ent.gz")
    t_path = os.path.join(pdb_dir, "pdb1a00.ent.gz")
    t_dat = os.path.join(DAT_DIR, "1a00A.dat")
    if not (os.path.exists(q_path) and os.path.exists(t_dat)):
        print(f"  align_pdb_with_dat: SKIPPED (files not found)")
        return

    # Align using .dat for template (bypasses Rust DSSP for template)
    result = dali.align_pdb(
        q_path, t_path, "A", "A", "101m", "1a00A",
        template_dat=t_dat,
    )
    assert result is not None, "expected AlignResult, got None"
    assert result.zscore > 0, f"expected zscore > 0, got {result.zscore}"
    assert result.n_aligned > 0, f"expected n_aligned > 0, got {result.n_aligned}"

    # Compare with PDB-only alignment
    result_pdb = dali.align_pdb(q_path, t_path, "A", "A", "101m", "1a00")
    assert result_pdb is not None

    # Both should find alignments (z-scores may differ due to different template DSSP)
    print(f"  align_pdb_with_dat: z={result.zscore:.1f} (dat), z={result_pdb.zscore:.1f} (pdb)")
    print(f"    n_aligned: {result.n_aligned} (dat), {result_pdb.n_aligned} (pdb)")


if __name__ == "__main__":
    tests = [
        ("read_dat", test_read_dat),
        ("ca_coords", test_ca_coords),
        ("import_pdb", test_import_pdb),
        ("ProteinStore", test_protein_store),
        ("compare_pair", test_compare_pair),
        ("run_wolf_path", test_run_wolf_path),
        ("resid_map", test_resid_map),
        ("write_dat_roundtrip", test_write_dat_roundtrip),
        ("align_pdb", test_align_pdb),
        ("align_pdb_with_dat", test_align_pdb_with_dat),
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
