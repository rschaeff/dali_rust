#!/usr/bin/env python3
"""Head-to-head benchmark: Rust dali-core vs Fortran DaliLite.v5.

Runs the same protein pairs through both backends and compares wall-clock times.
Methods: WOLF, PARSI, PIPELINE (Rust full), and the Fortran equivalents.

The PIPELINE method chains Fortran stages exactly as dali.pl does:
  Wolf path:  WOLF → DP → dccp2dalicon.pl → DALICON(lfitz=T) → DP
  Parsi path: PARSI → FILTER95 → sort|uniq|pipe96-free.pl → DALICON(lfitz=T) → DP
Both paths run in both directions (forward + reverse), matching Rust compare_pair.

Usage:
    python3 bench_compare.py [--methods WOLF,PARSI,PIPELINE] [--repeats 3]

Requires:
    - Rust binary: cargo build --release --features bench-bin --bin bench_head2head
    - Fortran: /home/rschaeff/src/Dali_v5/DaliLite.v5/bin/serialcompare + Perl scripts
"""

import argparse
import glob
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# --- Configuration ---

RUST_BINARY = Path(__file__).parent / "target" / "release" / "bench_head2head"
FORTRAN_BIN = Path("/home/rschaeff/src/Dali_v5/DaliLite.v5/bin")
FORTRAN_BINARY = FORTRAN_BIN / "serialcompare"
DCCP2DALICON = FORTRAN_BIN / "dccp2dalicon.pl"
PIPE96FREE = FORTRAN_BIN / "pipe96-free.pl"

DAT_DIR_5 = Path(__file__).parent.parent / "validation" / "fixtures" / "structures"
DAT_DIR_18 = Path(__file__).parent.parent / "validation" / "corpus_expansion" / "ground_truth_18" / "structures"

# Test pairs: (cd1, cd2, dat_dir, label, approx_combined_nres)
# Using 5-struct corpus only — Fortran segfaults on some 18-struct .dat files.
# All 5 structures: 101mA(154), 1a00A(141), 1a87A(297), 1allA(160), 1binA(143)
PAIRS = [
    ("101mA", "1a00A", DAT_DIR_5, "154+141", 295),
    ("101mA", "1binA", DAT_DIR_5, "154+143", 297),
    ("1a87A", "1allA", DAT_DIR_5, "297+160", 457),
    ("1a87A", "101mA", DAT_DIR_5, "297+154", 451),
    ("1a87A", "1binA", DAT_DIR_5, "297+143", 440),
    ("101mA", "1allA", DAT_DIR_5, "154+160", 314),
]


def ensure_trailing_slash(path):
    s = str(path)
    return s if s.endswith("/") else s + "/"


def run_rust(method, pairs, dat_dir, repeats=3):
    """Run Rust benchmark binary and return per-pair timings."""
    if not RUST_BINARY.exists():
        print(f"ERROR: Rust binary not found at {RUST_BINARY}")
        print("Build with: cd dali-core && cargo build --release --features bench-bin --bin bench_head2head")
        return None

    args = [str(RUST_BINARY), str(dat_dir), method]
    for cd1, cd2 in pairs:
        args.extend([cd1, cd2])

    best_results = None
    for r in range(repeats):
        result = subprocess.run(args, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  Rust error: {result.stderr.strip()}")
            return None
        data = json.loads(result.stdout)
        if best_results is None:
            best_results = data
        else:
            for i, pair in enumerate(data["pairs"]):
                if pair["elapsed_ms"] < best_results["pairs"][i]["elapsed_ms"]:
                    best_results["pairs"][i]["elapsed_ms"] = pair["elapsed_ms"]
            if data["total_ms"] < best_results["total_ms"]:
                best_results["total_ms"] = data["total_ms"]

    return best_results


def run_fortran_wolf(pairs, dat_dir, repeats=3):
    """Run Fortran serialcompare WOLF kernel only."""
    if not FORTRAN_BINARY.exists():
        print(f"ERROR: Fortran binary not found at {FORTRAN_BINARY}")
        return None

    dat_dir_str = ensure_trailing_slash(dat_dir)
    results = {"method": "WOLF", "backend": "fortran", "pairs": [], "total_ms": 0}

    total_start = time.time()
    for cd1, cd2 in pairs:
        best_ms = float("inf")
        n_results = 0
        for _ in range(repeats):
            with tempfile.TemporaryDirectory() as tmpdir:
                with open(os.path.join(tmpdir, "list1"), "w") as f:
                    f.write(f"{cd1}\n")
                with open(os.path.join(tmpdir, "list2"), "w") as f:
                    f.write(f"{cd2}\n")

                start = time.time()
                result = subprocess.run(
                    [str(FORTRAN_BINARY), dat_dir_str, dat_dir_str, "WOLF"],
                    cwd=tmpdir, capture_output=True, text=True, timeout=60,
                )
                elapsed_ms = (time.time() - start) * 1000

                if result.returncode != 0:
                    print(f"  Fortran WOLF error for {cd1}/{cd2}: {result.stderr.strip()[:200]}")
                    break

                fort101 = os.path.join(tmpdir, "fort.101")
                if os.path.exists(fort101):
                    with open(fort101) as f:
                        n_results = sum(1 for line in f if line.startswith("WOLFITZ"))

                best_ms = min(best_ms, elapsed_ms)

        results["pairs"].append({
            "cd1": cd1, "cd2": cd2,
            "elapsed_ms": round(best_ms, 3),
            "n_results": n_results,
        })

    results["total_ms"] = round((time.time() - total_start) * 1000 / repeats, 3)
    return results


def run_fortran_parsi(pairs, dat_dir, repeats=3):
    """Run Fortran serialcompare PARSI kernel only."""
    if not FORTRAN_BINARY.exists():
        return None

    dat_dir_str = ensure_trailing_slash(dat_dir)
    results = {"method": "PARSI", "backend": "fortran", "pairs": [], "total_ms": 0}
    total_start = time.time()

    for cd1, cd2 in pairs:
        best_ms = float("inf")
        n_results = 0
        for _ in range(repeats):
            with tempfile.TemporaryDirectory() as tmpdir:
                with open(os.path.join(tmpdir, "list1"), "w") as f:
                    f.write(f"{cd1}\n")
                with open(os.path.join(tmpdir, "list2"), "w") as f:
                    f.write(f"{cd2}\n")

                start = time.time()
                result = subprocess.run(
                    [str(FORTRAN_BINARY), dat_dir_str, dat_dir_str, "PARSI"],
                    cwd=tmpdir, capture_output=True, text=True, timeout=120,
                )
                elapsed_ms = (time.time() - start) * 1000

                if result.returncode != 0:
                    print(f"  Fortran PARSI error for {cd1}/{cd2}: {result.stderr.strip()[:200]}")
                    break

                fort101 = os.path.join(tmpdir, "fort.101")
                if os.path.exists(fort101):
                    with open(fort101) as f:
                        n_results = sum(1 for line in f if line.strip().startswith("refine"))

                best_ms = min(best_ms, elapsed_ms)

        results["pairs"].append({
            "cd1": cd1, "cd2": cd2,
            "elapsed_ms": round(best_ms, 3),
            "n_results": n_results,
        })

    results["total_ms"] = round((time.time() - total_start) * 1000 / repeats, 3)
    return results


def _collect_fort_output(tmpdir):
    """Collect all fort.1XX output files into a single string."""
    lines = []
    for f in sorted(glob.glob(os.path.join(tmpdir, "fort.1[0-9][0-9]"))):
        with open(f) as fh:
            lines.append(fh.read())
    return "".join(lines)


def _clean_fort_files(tmpdir):
    """Remove all fort.* files."""
    for f in glob.glob(os.path.join(tmpdir, "fort.*")):
        os.remove(f)


def _run_fortran_wolf_path(tmpdir, dat_dir_str):
    """Run the Fortran wolf path: WOLF → DP → dccp2dalicon.pl → DALICON(T) → DP.

    Mirrors mpidali.pm sub wolf() + sub dalicon().
    """
    sc = str(FORTRAN_BINARY)

    # Step 1: WOLF
    _clean_fort_files(tmpdir)
    r = subprocess.run([sc, dat_dir_str, dat_dir_str, "WOLF"],
                       cwd=tmpdir, capture_output=True, text=True, timeout=30)
    if r.returncode != 0:
        return False
    wolf_output = _collect_fort_output(tmpdir)
    if not wolf_output.strip():
        return True  # no WOLF hits — nothing to chain

    # Step 2: DP on WOLF output
    _clean_fort_files(tmpdir)
    r = subprocess.run([sc, dat_dir_str, dat_dir_str, "DP"],
                       input=wolf_output, cwd=tmpdir,
                       capture_output=True, text=True, timeout=30)
    if r.returncode != 0:
        return False
    dp_output = _collect_fort_output(tmpdir)
    if not dp_output.strip():
        return True

    # Step 3: dccp2dalicon.pl (Perl format converter)
    r = subprocess.run(["perl", str(DCCP2DALICON)],
                       input=dp_output, cwd=tmpdir,
                       capture_output=True, text=True, timeout=10)
    dalicon_input = r.stdout
    if not dalicon_input.strip():
        return True

    # Step 4: DALICON with lfitz=T
    _clean_fort_files(tmpdir)
    r = subprocess.run([sc, dat_dir_str, dat_dir_str, "DALICON", "T"],
                       input=dalicon_input, cwd=tmpdir,
                       capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        return False
    dalicon_output = _collect_fort_output(tmpdir)
    if not dalicon_output.strip():
        return True

    # Step 5: Final DP on DALICON output
    _clean_fort_files(tmpdir)
    r = subprocess.run([sc, dat_dir_str, dat_dir_str, "DP"],
                       input=dalicon_output, cwd=tmpdir,
                       capture_output=True, text=True, timeout=30)
    return r.returncode == 0


def _run_fortran_parsi_path(tmpdir, dat_dir_str):
    """Run the Fortran parsi path: PARSI → FILTER95 → pipe96-free.pl → DALICON(T) → DP.

    Mirrors mpidali.pm sub parsi() + sub dalicon().
    """
    sc = str(FORTRAN_BINARY)

    # Step 1: PARSI
    _clean_fort_files(tmpdir)
    r = subprocess.run([sc, dat_dir_str, dat_dir_str, "PARSI"],
                       cwd=tmpdir, capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        return False
    parsi_output = _collect_fort_output(tmpdir)
    if not parsi_output.strip():
        return True  # no PARSI hits

    # Step 2: FILTER95
    _clean_fort_files(tmpdir)
    r = subprocess.run([sc, dat_dir_str, dat_dir_str, "FILTER95"],
                       input=parsi_output, cwd=tmpdir,
                       capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        return False
    filter_output = _collect_fort_output(tmpdir)
    if not filter_output.strip():
        return True

    # Step 3: sort -nr | uniq | pipe96-free.pl 1.0 2
    r = subprocess.run(
        f"sort -nr | uniq | perl {PIPE96FREE} 1.0 2",
        input=filter_output, cwd=tmpdir, shell=True,
        capture_output=True, text=True, timeout=10,
    )
    dalicon_input = r.stdout
    if not dalicon_input.strip():
        return True

    # Step 4: DALICON with lfitz=T
    _clean_fort_files(tmpdir)
    r = subprocess.run([sc, dat_dir_str, dat_dir_str, "DALICON", "T"],
                       input=dalicon_input, cwd=tmpdir,
                       capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        return False
    dalicon_output = _collect_fort_output(tmpdir)
    if not dalicon_output.strip():
        return True

    # Step 5: Final DP on DALICON output
    _clean_fort_files(tmpdir)
    r = subprocess.run([sc, dat_dir_str, dat_dir_str, "DP"],
                       input=dalicon_output, cwd=tmpdir,
                       capture_output=True, text=True, timeout=30)
    return r.returncode == 0


def run_fortran_pipeline(pairs, dat_dir, repeats=3):
    """Run the full Fortran pipeline, chained exactly as dali.pl does.

    For each pair, runs:
      Forward:  wolf_path(cd1→cd2) + parsi_path(cd1→cd2)
      Reverse:  wolf_path(cd2→cd1) + parsi_path(cd2→cd1)

    This matches what Rust compare_pair() does.
    """
    if not FORTRAN_BINARY.exists():
        print(f"ERROR: Fortran binary not found at {FORTRAN_BINARY}")
        return None

    dat_dir_str = ensure_trailing_slash(dat_dir)
    results = {"method": "PIPELINE", "backend": "fortran", "pairs": [], "total_ms": 0}
    total_start = time.time()

    for cd1, cd2 in pairs:
        best_ms = float("inf")
        for _ in range(repeats):
            with tempfile.TemporaryDirectory() as tmpdir:
                start = time.time()

                # Forward direction: cd1 → cd2
                with open(os.path.join(tmpdir, "list1"), "w") as f:
                    f.write(f"{cd1}\n")
                with open(os.path.join(tmpdir, "list2"), "w") as f:
                    f.write(f"{cd2}\n")
                _run_fortran_wolf_path(tmpdir, dat_dir_str)
                _run_fortran_parsi_path(tmpdir, dat_dir_str)

                # Reverse direction: cd2 → cd1
                with open(os.path.join(tmpdir, "list1"), "w") as f:
                    f.write(f"{cd2}\n")
                with open(os.path.join(tmpdir, "list2"), "w") as f:
                    f.write(f"{cd1}\n")
                _run_fortran_wolf_path(tmpdir, dat_dir_str)
                _run_fortran_parsi_path(tmpdir, dat_dir_str)

                elapsed_ms = (time.time() - start) * 1000
                best_ms = min(best_ms, elapsed_ms)

        results["pairs"].append({
            "cd1": cd1, "cd2": cd2,
            "elapsed_ms": round(best_ms, 3),
            "n_results": 0,
        })

    results["total_ms"] = round((time.time() - total_start) * 1000 / repeats, 3)
    return results


def print_comparison(rust_data, fortran_data, method):
    """Print a side-by-side comparison table."""
    if rust_data is None or fortran_data is None:
        print(f"\n  [{method}] Skipped — one backend failed\n")
        return

    print(f"\n{'='*80}")
    print(f"  {method} — Rust vs Fortran (best of N)")
    print(f"{'='*80}")
    print(f"  {'Pair':<24} {'Nres':>6} {'Rust (ms)':>12} {'Fortran (ms)':>14} {'Ratio':>8} {'Winner':>8}")
    print(f"  {'-'*24} {'-'*6} {'-'*12} {'-'*14} {'-'*8} {'-'*8}")

    for i, (rp, fp) in enumerate(zip(rust_data["pairs"], fortran_data["pairs"])):
        label = f"{rp['cd1']}/{rp['cd2']}"
        nres = PAIRS[i][4] if i < len(PAIRS) else 0
        rust_ms = rp["elapsed_ms"]
        fort_ms = fp["elapsed_ms"]
        ratio = rust_ms / fort_ms if fort_ms > 0 else float("inf")
        winner = "Rust" if ratio < 1.0 else "Fortran"
        print(f"  {label:<24} {nres:>6} {rust_ms:>12.1f} {fort_ms:>14.1f} {ratio:>8.2f}x {winner:>8}")

    rust_total = sum(p["elapsed_ms"] for p in rust_data["pairs"])
    fort_total = sum(p["elapsed_ms"] for p in fortran_data["pairs"])
    ratio = rust_total / fort_total if fort_total > 0 else float("inf")
    print(f"  {'-'*24} {'-'*6} {'-'*12} {'-'*14} {'-'*8} {'-'*8}")
    print(f"  {'TOTAL':<24} {'':>6} {rust_total:>12.1f} {fort_total:>14.1f} {ratio:>8.2f}x")
    print()


def main():
    parser = argparse.ArgumentParser(description="Rust vs Fortran DALI benchmark")
    parser.add_argument("--methods", default="WOLF,PARSI,PIPELINE",
                        help="Comma-separated methods to benchmark")
    parser.add_argument("--repeats", type=int, default=3,
                        help="Number of repeats (best-of-N)")
    parser.add_argument("--json", default=None,
                        help="Save raw results to JSON file")
    args = parser.parse_args()

    methods = [m.strip() for m in args.methods.split(",")]

    print(f"Benchmark: Rust dali-core vs Fortran DaliLite.v5")
    print(f"Repeats: {args.repeats} (best-of-N)")
    print(f"Methods: {', '.join(methods)}")
    print(f"Rust binary: {RUST_BINARY}")
    print(f"Fortran binary: {FORTRAN_BINARY}")

    all_results = {}

    for method in methods:
        print(f"\nRunning {method}...")

        # Group pairs by dat_dir
        groups = {}
        for cd1, cd2, dat_dir, label, nres in PAIRS:
            key = str(dat_dir)
            if key not in groups:
                groups[key] = []
            groups[key].append((cd1, cd2))

        # Rust
        rust_combined = {"method": method, "backend": "rust", "pairs": [], "total_ms": 0}
        for dat_dir_str, group_pairs in groups.items():
            r = run_rust(method, group_pairs, dat_dir_str, args.repeats)
            if r:
                rust_combined["pairs"].extend(r["pairs"])
                rust_combined["total_ms"] += r["total_ms"]

        # Fortran
        fortran_data = None
        if method == "WOLF":
            fort_combined = {"method": "WOLF", "backend": "fortran", "pairs": [], "total_ms": 0}
            for dat_dir_str, group_pairs in groups.items():
                f = run_fortran_wolf(group_pairs, dat_dir_str, args.repeats)
                if f:
                    fort_combined["pairs"].extend(f["pairs"])
                    fort_combined["total_ms"] += f["total_ms"]
            fortran_data = fort_combined
        elif method == "PARSI":
            fort_combined = {"method": "PARSI", "backend": "fortran", "pairs": [], "total_ms": 0}
            for dat_dir_str, group_pairs in groups.items():
                f = run_fortran_parsi(group_pairs, dat_dir_str, args.repeats)
                if f:
                    fort_combined["pairs"].extend(f["pairs"])
                    fort_combined["total_ms"] += f["total_ms"]
            fortran_data = fort_combined
        elif method == "PIPELINE":
            fort_combined = {"method": "PIPELINE", "backend": "fortran", "pairs": [], "total_ms": 0}
            for dat_dir_str, group_pairs in groups.items():
                f = run_fortran_pipeline(group_pairs, dat_dir_str, args.repeats)
                if f:
                    fort_combined["pairs"].extend(f["pairs"])
                    fort_combined["total_ms"] += f["total_ms"]
            fortran_data = fort_combined

        rust_data = rust_combined if rust_combined["pairs"] else None
        print_comparison(rust_data, fortran_data, method)

        all_results[method] = {"rust": rust_data, "fortran": fortran_data}

    if args.json:
        with open(args.json, "w") as f:
            json.dump(all_results, f, indent=2)
        print(f"Raw results saved to {args.json}")


if __name__ == "__main__":
    main()
