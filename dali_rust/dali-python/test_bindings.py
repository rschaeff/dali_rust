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


if __name__ == "__main__":
    tests = [
        ("read_dat", test_read_dat),
        ("ca_coords", test_ca_coords),
        ("import_pdb", test_import_pdb),
        ("ProteinStore", test_protein_store),
        ("compare_pair", test_compare_pair),
        ("run_wolf_path", test_run_wolf_path),
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
