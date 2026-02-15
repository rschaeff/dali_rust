#!/usr/bin/env python3
"""Test FILTER95 module against 5-structure corpus ground truth."""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

# Find .dat files
STRUCTURE_DAT_DIR = None
for candidate in [
    Path('/home/rschaeff/src/Dali_v5/DaliLite.v5/DAT'),
    Path('validation/fixtures/wolf/structures'),
]:
    if candidate.exists():
        STRUCTURE_DAT_DIR = candidate
        break

if STRUCTURE_DAT_DIR is None:
    print("ERROR: Could not find .dat files directory")
    sys.exit(1)

print(f"Using .dat files from: {STRUCTURE_DAT_DIR}")

GROUND_TRUTH_INPUT = Path('validation/fixtures/filter95/filter95_input.txt')
GROUND_TRUTH_OUTPUT = Path('validation/fixtures/filter95/filter95_output.txt')


def test_smoke():
    """Quick smoke test: parse one refine line and process it."""
    from filter95.filter95 import parse_refine_line, readproteindata95

    print("=" * 70)
    print("FILTER95 Module - Smoke test")
    print("=" * 70)

    # Test refine line parsing
    with open(GROUND_TRUTH_INPUT) as f:
        first_line = f.readline()

    parsed = parse_refine_line(first_line)
    if parsed is None:
        print("  FAILED to parse refine line")
        return False

    cd1, cd2, idom, score, nseg, ali = parsed
    print(f"  Parsed: cd1={cd1} cd2={cd2} idom={idom} score={score} nseg={nseg}")
    print(f"  Ali ({len(ali)} values): {ali[:8]}...")

    # Test .dat reader with domain info
    dat_path = str(STRUCTURE_DAT_DIR) + '/101mA.dat'
    nres, ca, ndom, node_type, node_size, node_nseg, segrange = \
        readproteindata95(dat_path)
    print(f"  101mA: nres={nres} ndom={ndom}")
    for i in range(1, ndom + 1):
        print(f"    dom {i}: type={node_type[i]} size={node_size[i]} "
              f"nseg={node_nseg[i]} ranges={segrange.get(i, [])}")

    return True


def test_single_pair():
    """Test FILTER95 on a small subset of refine lines."""
    from filter95.filter95 import run_filter95

    print("\n" + "=" * 70)
    print("FILTER95 Module - Single pair test (101mA vs 101mA)")
    print("=" * 70)

    # Get just the 101mA101mA refine lines
    with open(GROUND_TRUTH_INPUT) as f:
        all_lines = f.readlines()

    pair_lines = [l for l in all_lines if '101mA101mA' in l]
    print(f"  Input: {len(pair_lines)} refine lines for 101mA101mA")

    t0 = time.time()
    result = run_filter95(pair_lines, str(STRUCTURE_DAT_DIR) + '/')
    elapsed = time.time() - t0

    print(f"  Output: {len(result)} lines in {elapsed:.2f}s")
    for line in result[:3]:
        # Show first ~80 chars
        print(f"    {line[:80]}...")

    return len(result) > 0


def test_full_corpus():
    """Run FILTER95 on full 5-structure corpus and compare to ground truth."""
    from filter95.filter95 import run_filter95

    print("\n" + "=" * 70)
    print("FILTER95 Module - Full corpus (5 structures)")
    print("=" * 70)

    # Load input (PARSI refine lines)
    with open(GROUND_TRUTH_INPUT) as f:
        input_lines = f.readlines()
    print(f"  Input: {len(input_lines)} refine lines")

    # Load ground truth output
    with open(GROUND_TRUTH_OUTPUT) as f:
        ref_text = f.read()
    ref_lines = [l.strip() for l in ref_text.strip().split('\n') if l.strip()]
    print(f"  Reference: {len(ref_lines)} output lines")

    # Run FILTER95
    t0 = time.time()
    result_lines = run_filter95(input_lines, str(STRUCTURE_DAT_DIR) + '/')
    elapsed = time.time() - t0
    print(f"  Python: {len(result_lines)} output lines ({elapsed:.1f}s)")

    # Parse and compare
    ref_parsed = _parse_output(ref_lines)
    cand_parsed = _parse_output(result_lines)

    matched = 0
    mismatched = 0

    for key, ref in ref_parsed.items():
        cand = cand_parsed.get(key)
        if cand is None:
            mismatched += 1
            if mismatched <= 10:
                print(f"  MISSING: {key} zscore={ref['zscore']:.3f} "
                      f"score={ref['score']} nres={ref['nres']}")
            continue

        # Compare: zscore within tolerance, score within tolerance
        z_ok = abs(ref['zscore'] - cand['zscore']) < 0.5
        s_ok = abs(ref['score'] - cand['score']) <= max(abs(ref['score']) * 0.01, 100)
        n_ok = abs(ref['nres'] - cand['nres']) <= 5

        if z_ok and s_ok:
            matched += 1
        else:
            mismatched += 1
            if mismatched <= 10:
                print(f"  MISMATCH: {key}")
                print(f"    ref:  z={ref['zscore']:.3f} score={ref['score']} "
                      f"nres={ref['nres']}")
                print(f"    cand: z={cand['zscore']:.3f} score={cand['score']} "
                      f"nres={cand['nres']}")

    # Check for extra lines
    extra = len(cand_parsed) - len(ref_parsed)

    print(f"\n  FILTER95: {matched}/{len(ref_parsed)} matched, "
          f"{mismatched} missing/mismatched, {extra} extra")

    return matched, mismatched


def _parse_output(lines):
    """Parse FILTER95 output lines into dict keyed by (cd1cd2, idom)."""
    parsed = {}
    for line in lines:
        parts = line.split()
        if len(parts) < 5:
            continue
        zscore = float(parts[0])
        cd1cd2 = parts[1]
        idom = int(parts[2])
        score = int(parts[3])
        nres = int(parts[4])

        key = (cd1cd2, idom)
        parsed[key] = {
            'zscore': zscore,
            'score': score,
            'nres': nres,
        }
    return parsed


if __name__ == '__main__':
    if '--smoke' in sys.argv:
        ok = test_smoke()
        sys.exit(0 if ok else 1)
    elif '--single' in sys.argv:
        ok = test_single_pair()
        sys.exit(0 if ok else 1)
    elif '--full' in sys.argv:
        matched, mismatched = test_full_corpus()
        sys.exit(0 if mismatched == 0 else 1)
    else:
        ok = test_smoke()
        if ok:
            ok = test_single_pair()
        if ok:
            test_full_corpus()
