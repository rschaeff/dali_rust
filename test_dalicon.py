#!/usr/bin/env python3
"""
Test DALICON module against ground truth fixtures.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from dalicon.dalicon import run_dalicon
from validation.parsers import parse_wolfitz


DAT_DIR = Path('validation/fixtures/structures')
DALICON_INPUT = Path('validation/fixtures/dalicon/wolf/dalicon_input.txt')
DALICON_FIXTURE = Path('validation/fixtures/dalicon/wolf/dalicon_output.txt')


def main():
    print("=" * 60)
    print("DALICON Module Validation Test")
    print("=" * 60)

    # Step 1: Run DALICON
    print("\n--- Step 1: Running DALICON ---")
    output_lines = run_dalicon(DALICON_INPUT, DAT_DIR)
    print(f"  Produced {len(output_lines)} WOLFITZ hits")

    # Write output
    output_path = Path('dalicon_output_new.txt')
    with open(output_path, 'w') as f:
        for line in output_lines:
            f.write(line + '\n')
    print(f"  Written to {output_path}")

    # Step 2: Compare against ground truth
    print("\n--- Step 2: Comparison against ground truth ---")
    ref_hits = parse_wolfitz(DALICON_FIXTURE)
    ref_map = {}
    for h in ref_hits:
        key = (h.cd1, h.cd2)
        ref_map[key] = h

    # Parse our output
    cand_hits = parse_wolfitz(output_path)
    cand_map = {}
    for h in cand_hits:
        key = (h.cd1, h.cd2)
        cand_map[key] = h

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
            print(f"  EXTRA:   {pair_label}")
            extra += 1
        else:
            ref = ref_map[key]
            cand = cand_map[key]

            issues = []
            if ref.nblock != cand.nblock:
                issues.append(f"nblock: ref={ref.nblock} cand={cand.nblock}")
            if ref.raw_values != cand.raw_values:
                # Show first difference
                for idx, (rv, cv) in enumerate(zip(ref.raw_values, cand.raw_values)):
                    if rv != cv:
                        issues.append(f"value[{idx}]: ref={rv} cand={cv}")
                        break
                if len(ref.raw_values) != len(cand.raw_values):
                    issues.append(f"value count: ref={len(ref.raw_values)} cand={len(cand.raw_values)}")

            if issues:
                print(f"  MISMATCH: {pair_label}: {'; '.join(issues)}")
                mismatched += 1
            else:
                print(f"  MATCH:   {pair_label} (nblock={cand.nblock})")
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
