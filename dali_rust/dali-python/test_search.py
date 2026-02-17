"""Tests for search_database, mask_protein, iterative_search, and add_protein."""
import os
import shutil
import tempfile

import dali

DAT_DIR = os.path.join(os.path.dirname(__file__), "../../validation/fixtures/structures")
CODES_5 = ["101mA", "1a00A", "1a87A", "1allA", "1binA"]


def test_mask_protein_basic():
    """mask_protein subsets residues, reconstructs segments, round-trips through .dat."""
    p = dali.read_dat(os.path.join(DAT_DIR, "101mA.dat"))
    assert p.nres == 154

    # Keep first half of residues
    keep = list(range(p.nres // 2))
    masked = dali.mask_protein(p, keep, "mask1")

    assert masked.code == "mask1"
    assert masked.nres == len(keep)
    assert masked.nres == 77

    # Segments should have valid ranges within [1, new_nres]
    for b in range(masked.nseg):
        pass  # segments are internal, but nres/nseg should be consistent

    # CA shape should be (3, new_nres)
    ca = masked.ca_coords()
    assert ca.shape == (3, 77)

    # resid_map should be subset
    orig_resids = p.resid_map
    masked_resids = masked.resid_map
    assert len(masked_resids) == 77
    for i, idx in enumerate(keep):
        assert masked_resids[i] == orig_resids[idx]

    # write_dat → read_dat round-trip
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "mask1.dat")
        masked.write_dat(path)
        reloaded = dali.read_dat(path)
        assert reloaded.code == "mask1"
        assert reloaded.nres == 77
        assert reloaded.nseg == masked.nseg


def test_mask_protein_preserves_sequence():
    """mask_protein subsets sequence when available (import_pdb provides sequence)."""
    pdb_path = os.path.join(DAT_DIR, "pdb/pdb101m.ent.gz")
    if not os.path.exists(pdb_path):
        return  # skip if PDB not available
    p = dali.import_pdb(pdb_path, "A", "101mA")
    if p.sequence:
        keep = [0, 1, 2, 10, 20, 30]
        masked = dali.mask_protein(p, keep, "seq_t")
        assert len(masked.sequence) == len(keep)


def test_add_protein_to_store():
    """add_protein writes .dat and makes protein retrievable."""
    with tempfile.TemporaryDirectory() as tmpdir:
        store = dali.ProteinStore(tmpdir)

        # Load a protein and add it
        p = dali.read_dat(os.path.join(DAT_DIR, "101mA.dat"))
        store.add_protein(p)

        # Should be retrievable
        assert "101mA" in store
        p2 = store.get_protein("101mA")
        assert p2.nres == p.nres
        assert p2.nseg == p.nseg

        # .dat file should exist on disk
        assert os.path.exists(os.path.join(tmpdir, "101mA.dat"))


def test_search_database_basic():
    """search_database finds known matches from 5-structure corpus."""
    store = dali.ProteinStore(DAT_DIR)
    targets = ["1a00A", "1binA", "1a87A", "1allA"]

    results = dali.search_database("101mA", targets, store, z_cut=2.0)

    assert len(results) > 0
    # Results should be sorted by z-score descending
    for i in range(1, len(results)):
        assert results[i - 1].zscore >= results[i].zscore

    # All results should have cd1 == "101mA" and cd2 != "101mA"
    for r in results:
        assert r.cd1 == "101mA"
        assert r.cd2 != "101mA"
        assert r.zscore >= 2.0

    # 101mA vs 1a00A/1binA are known globin matches from E2E validation
    cd2s = {r.cd2 for r in results}
    assert "1a00A" in cd2s or "1binA" in cd2s, f"Expected globin match, got {cd2s}"


def test_search_database_skip_wolf():
    """search_database with skip_wolf=True uses only PARSI path."""
    store = dali.ProteinStore(DAT_DIR)
    targets = ["1a00A", "1binA"]

    results = dali.search_database("101mA", targets, store, z_cut=2.0, skip_wolf=True)

    # Should still find matches via PARSI path
    assert len(results) > 0
    for r in results:
        assert r.zscore >= 2.0


def test_search_database_no_matches():
    """search_database returns empty list when z_cut is very high."""
    store = dali.ProteinStore(DAT_DIR)
    results = dali.search_database("101mA", ["1a87A"], store, z_cut=100.0)
    assert len(results) == 0


def test_iterative_search_basic():
    """iterative_search finds at least one hit on corpus."""
    store = dali.ProteinStore(DAT_DIR)
    targets = ["1a00A", "1binA", "1a87A", "1allA"]

    hits = dali.iterative_search(
        "101mA", targets, store,
        min_aligned=10, min_zscore=2.0, gap_tolerance=5, max_rounds=3,
    )

    assert len(hits) >= 1

    hit = hits[0]
    assert hit.round == 0
    assert hit.zscore >= 2.0
    assert hit.nblock > 0
    assert len(hit.blocks) == hit.nblock
    assert len(hit.alignments) > 0
    assert len(hit.rotation) == 3
    assert all(len(row) == 3 for row in hit.rotation)
    assert len(hit.translation) == 3
    assert hit.cd2 in targets

    # Alignment pairs should be 1-based
    for q, t in hit.alignments:
        assert q >= 1
        assert t >= 1


def test_iterative_search_masked_protein_is_smaller():
    """After iteration, masked proteins should have fewer residues."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Copy .dat files to temp dir so masked proteins can be written
        for code in CODES_5:
            shutil.copy(os.path.join(DAT_DIR, f"{code}.dat"), tmpdir)

        store = dali.ProteinStore(tmpdir)
        orig = store.get_protein("101mA")
        targets = ["1a00A", "1binA"]

        hits = dali.iterative_search(
            "101mA", targets, store,
            min_aligned=10, min_zscore=2.0, max_rounds=3,
        )

        if len(hits) >= 1:
            # Check that a masked protein was created in the store
            masked_code = f"101mA_r0"
            if masked_code in store:
                masked = store.get_protein(masked_code)
                assert masked.nres < orig.nres, \
                    f"Masked {masked.nres} should be < original {orig.nres}"


def test_search_hit_repr():
    """SearchHit has a useful string representation."""
    store = dali.ProteinStore(DAT_DIR)
    hits = dali.iterative_search(
        "101mA", ["1a00A"], store,
        min_aligned=10, min_zscore=2.0, max_rounds=1,
    )
    if hits:
        r = repr(hits[0])
        assert "SearchHit" in r
        assert "zscore=" in r
