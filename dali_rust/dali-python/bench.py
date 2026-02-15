"""Benchmark script for dali PyO3 bindings. Run after each optimization."""

import time
import sys
import dali

import os
_HERE = os.path.dirname(os.path.abspath(__file__))
DAT_DIR_5 = os.path.join(_HERE, "../../validation/fixtures/structures")
DAT_DIR_18 = os.path.join(_HERE, "../../validation/corpus_expansion/ground_truth_18/structures")

CODES_5 = ["101mA", "1a00A", "1a87A", "1allA", "1binA"]


def bench_single_pair():
    """Time a single compare_pair call (101mA vs 1a00A)."""
    store = dali.ProteinStore(DAT_DIR_5)
    # Warm the cache
    store.get_protein("101mA")
    store.get_protein("1a00A")

    times = []
    for _ in range(3):
        t0 = time.perf_counter()
        results = dali.compare_pair("101mA", "1a00A", store)
        t1 = time.perf_counter()
        times.append(t1 - t0)
    best = min(times)
    print(f"single_pair:  {best:.3f}s  ({len(results)} results, best of 3)")
    return best


def bench_20_pairs():
    """Time all 20 directed pairs from 5 structures."""
    store = dali.ProteinStore(DAT_DIR_5)
    # Warm the cache
    for c in CODES_5:
        store.get_protein(c)

    t0 = time.perf_counter()
    total = 0
    for c1 in CODES_5:
        for c2 in CODES_5:
            if c1 != c2:
                results = dali.compare_pair(c1, c2, store)
                total += len(results)
    t1 = time.perf_counter()
    elapsed = t1 - t0
    print(f"20_pairs:     {elapsed:.3f}s  ({total} results, {elapsed/20:.3f}s/pair)")
    return elapsed


def bench_wolf_only():
    """Time just the WOLF path for all 20 pairs."""
    store = dali.ProteinStore(DAT_DIR_5)
    for c in CODES_5:
        store.get_protein(c)

    t0 = time.perf_counter()
    total = 0
    for c1 in CODES_5:
        targets = [c2 for c2 in CODES_5 if c2 != c1]
        results = dali.run_wolf_path(c1, targets, store)
        total += len(results)
    t1 = time.perf_counter()
    elapsed = t1 - t0
    print(f"wolf_20:      {elapsed:.3f}s  ({total} results)")
    return elapsed


if __name__ == "__main__":
    print("=== dali benchmark ===")
    bench_single_pair()
    bench_20_pairs()
    bench_wolf_only()
