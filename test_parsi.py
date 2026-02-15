#!/usr/bin/env python3
"""Test PARSI module against 5-structure corpus ground truth."""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

DAT_DIR = Path('validation/fixtures/parsi/')
# Try to find .dat files — they should be in the same place as other fixtures
STRUCTURE_DAT_DIR = None
for candidate in [
    Path('/home/rschaeff/src/Dali_v5/DaliLite.v5/DAT'),
    Path('validation/fixtures/wolf/structures'),
    Path('validation/corpus_expansion/ground_truth_18/structures'),
]:
    if candidate.exists():
        STRUCTURE_DAT_DIR = candidate
        break

if STRUCTURE_DAT_DIR is None:
    print("ERROR: Could not find .dat files directory")
    sys.exit(1)

print(f"Using .dat files from: {STRUCTURE_DAT_DIR}")

STRUCTURES = ['101mA', '1a00A', '1a87A', '1allA', '1binA']
GROUND_TRUTH = Path('validation/fixtures/parsi/parsi_output.txt')


def test_single_pair():
    """Quick smoke test: just run one pair."""
    from parsi.parsi import dowork_parsi, init_parsi

    print("=" * 70)
    print("PARSI Module - Single pair smoke test (101mA vs 101mA)")
    print("=" * 70)

    init_parsi()
    t0 = time.time()
    lines = dowork_parsi('101mA', '101mA', str(STRUCTURE_DAT_DIR) + '/')
    elapsed = time.time() - t0

    print(f"  Got {len(lines)} refine lines in {elapsed:.1f}s")
    for line in lines[:5]:
        print(f"    {line.strip()}")
    if len(lines) > 5:
        print(f"    ... ({len(lines) - 5} more)")

    return len(lines) > 0


def test_full_corpus():
    """Run PARSI on all 5x5 pairs and compare to ground truth."""
    from parsi.parsi import run_parsi

    print("\n" + "=" * 70)
    print("PARSI Module - Full corpus (5 structures, 25 pairs)")
    print("=" * 70)

    # Load ground truth
    ref_lines = []
    with open(GROUND_TRUTH) as f:
        for line in f:
            line = line.strip()
            if line:
                ref_lines.append(line)
    print(f"  Reference: {len(ref_lines)} refine lines")

    # Run PARSI
    t0 = time.time()
    result_lines = run_parsi(STRUCTURES, str(STRUCTURE_DAT_DIR) + '/')
    elapsed = time.time() - t0
    print(f"  Python: {len(result_lines)} refine lines ({elapsed:.1f}s)")

    # Parse and compare
    ref_parsed = _parse_refine_lines(ref_lines)
    cand_parsed = _parse_refine_lines(result_lines)

    ref_map = {}
    for r in ref_parsed:
        key = (r['cd1cd2'], r['idom'])
        ref_map.setdefault(key, []).append(r)

    cand_map = {}
    for c in cand_parsed:
        key = (c['cd1cd2'], c['idom'])
        cand_map.setdefault(key, []).append(c)

    matched = 0
    mismatched = 0
    missing = 0

    for key in sorted(ref_map.keys()):
        ref_group = ref_map[key]
        cand_group = cand_map.get(key, [])

        for ref in ref_group:
            found = False
            for cand in cand_group:
                if (ref['score'] == cand['score'] and
                        ref['nseg'] == cand['nseg'] and
                        ref['ranges'] == cand['ranges']):
                    found = True
                    break
            if found:
                matched += 1
            else:
                # Check with score tolerance
                for cand in cand_group:
                    if (abs(ref['score'] - cand['score']) <= 100 and
                            ref['nseg'] == cand['nseg']):
                        found = True
                        break
                if found:
                    matched += 1
                else:
                    mismatched += 1
                    if mismatched <= 10:
                        print(f"  MISSING: {key} score={ref['score']} nseg={ref['nseg']}")

    print(f"\n  PARSI: {matched}/{matched+mismatched} matched, "
          f"{mismatched} missing")

    return matched, mismatched


def _parse_refine_lines(lines):
    """Parse refine output lines."""
    parsed = []
    for line in lines:
        line = line.strip()
        if not line or 'refine' not in line:
            continue
        parts = line.split()
        # First part: refine<cd1cd2>
        cd1cd2 = parts[0].replace('refine', '')
        idom = int(parts[1])
        score = int(parts[2])
        nseg = int(parts[3])
        ranges = [int(x) for x in parts[4:]]
        parsed.append({
            'cd1cd2': cd1cd2,
            'idom': idom,
            'score': score,
            'nseg': nseg,
            'ranges': ranges,
        })
    return parsed


if __name__ == '__main__':
    if '--smoke' in sys.argv:
        ok = test_single_pair()
        sys.exit(0 if ok else 1)
    elif '--full' in sys.argv:
        matched, mismatched = test_full_corpus()
        sys.exit(0 if mismatched == 0 else 1)
    else:
        # Default: smoke test first
        ok = test_single_pair()
        if ok:
            test_full_corpus()
