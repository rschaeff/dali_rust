#!/usr/bin/env python3
"""
Test WOLF module against ground truth fixtures.
"""

import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))

from wolf.wolf import run_wolf_all_pairs, WolfResult
from validation.parsers import parse_wolfitz


DAT_DIR = Path('validation/fixtures/structures')
WOLF_FIXTURE = Path('validation/fixtures/wolf/wolf_output.txt')

STRUCTURES = ['101mA', '1a00A', '1a87A', '1allA', '1binA']


def main():
    print("=" * 60)
    print("WOLF Module Validation Test")
    print("=" * 60)

    # Step 1: Basic .dat reading test
    print("\n--- Step 1: .dat file reading ---")
    from wolf.dat_reader import read_dat
    for code in STRUCTURES:
        p = read_dat(DAT_DIR / f'{code}.dat')
        print(f"  {p.code}: nres={p.nres}, nseg={p.nseg}, "
              f"na={p.na}, nb={p.nb}, secstr={''.join(p.secstr)}")

    # Step 2: Run WOLF on all pairs
    print("\n--- Step 2: Running WOLF comparisons ---")
    results = run_wolf_all_pairs(STRUCTURES, DAT_DIR)
    print(f"  Produced {len(results)} WOLFITZ hits")

    # Step 3: Write output
    output_path = Path('wolf_output_new.txt')
    with open(output_path, 'w') as f:
        for r in results:
            f.write(r.to_wolfitz_line() + '\n')
    print(f"  Written to {output_path}")

    # Step 4: Compare against ground truth
    print("\n--- Step 3: Comparison against ground truth ---")
    ref_hits = parse_wolfitz(WOLF_FIXTURE)
    ref_map = {(h.cd1, h.cd2): h for h in ref_hits}

    cand_map = {(r.cd1, r.cd2): r for r in results}

    all_keys = sorted(set(ref_map.keys()) | set(cand_map.keys()))
    matched = 0
    mismatched = 0
    missing = 0
    extra = 0

    for key in all_keys:
        cd1, cd2 = key
        pair_label = f'{cd1}-{cd2}'

        if key not in cand_map:
            print(f"  MISSING: {pair_label}")
            missing += 1
        elif key not in ref_map:
            print(f"  EXTRA:   {pair_label} (nblock={cand_map[key].nblock})")
            extra += 1
        else:
            ref_h = ref_map[key]
            cand_r = cand_map[key]

            # Build raw_values from our result for comparison
            cand_vals = []
            for i in range(cand_r.nblock):
                cand_vals.append(cand_r.l1[i])
                cand_vals.append(cand_r.r1[i])
            for i in range(cand_r.nblock):
                cand_vals.append(cand_r.l2[i])
                cand_vals.append(cand_r.r2[i])

            if ref_h.nblock != cand_r.nblock:
                print(f"  MISMATCH: {pair_label} nblock: ref={ref_h.nblock} cand={cand_r.nblock}")
                mismatched += 1
            elif ref_h.raw_values != cand_vals:
                # Show first difference
                for i, (rv, cv) in enumerate(zip(ref_h.raw_values, cand_vals)):
                    if rv != cv:
                        print(f"  MISMATCH: {pair_label} value[{i}]: ref={rv} cand={cv}")
                        break
                mismatched += 1
            else:
                print(f"  MATCH:   {pair_label} (nblock={cand_r.nblock})")
                matched += 1

    print(f"\n{'=' * 60}")
    print(f"Results: {matched} matched, {mismatched} mismatched, "
          f"{missing} missing, {extra} extra")
    total = matched + mismatched + missing
    if mismatched == 0 and missing == 0:
        print("ALL PASSED!")
    else:
        print(f"FAILED: {mismatched + missing} issues")
    print(f"{'=' * 60}")

    return 0 if (mismatched == 0 and missing == 0) else 1


if __name__ == '__main__':
    sys.exit(main())
