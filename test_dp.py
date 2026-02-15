#!/usr/bin/env python3
"""
Test DP module against ground truth fixtures.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from dp.dp import run_dp, dpsetup
from validation.parsers import parse_dccp


DAT_DIR = Path('validation/fixtures/structures')
DP_INPUT = Path('validation/fixtures/dp/wolf/dp_input.txt')
DP_FIXTURE = Path('validation/fixtures/dp/wolf/dp_output.txt')


def main():
    print("=" * 60)
    print("DP Module Validation Test")
    print("=" * 60)

    # Step 1: Test dpsetup
    print("\n--- Step 1: dpsetup ---")
    for code in ['101mA', '1a00A', '1a87A', '1allA', '1binA']:
        nres, ca, d, ndom, domns, domseglist = dpsetup(DAT_DIR / f'{code}.dat')
        print(f"  {code}: nres={nres}, ndom={ndom}")
        for i in range(ndom):
            segs_str = ', '.join(f'{s[0]}-{s[1]}' for s in domseglist[i])
            print(f"    domain {i}: {domns[i]} seg(s) [{segs_str}]")

    # Step 2: Run DP
    print("\n--- Step 2: Running DP ---")
    results = run_dp(DP_INPUT, DAT_DIR)
    print(f"  Produced {len(results)} DCCP hits")

    # Step 3: Write output
    output_path = Path('dp_output_new.txt')
    with open(output_path, 'w') as f:
        for r in results:
            f.write(r.format_dccp() + '\n')
    print(f"  Written to {output_path}")

    # Step 4: Compare against ground truth
    print("\n--- Step 3: Comparison against ground truth ---")
    ref_blocks = parse_dccp(DP_FIXTURE)
    ref_map = {}
    for b in ref_blocks:
        key = (b.cd1, b.cd2)
        ref_map[key] = b

    cand_map = {}
    for r in results:
        key = (r.cd1, r.cd2)
        # Build DCCPBlock-compatible ranges
        ranges_cd1 = [(r.l1[i], r.r1[i]) for i in range(r.nblock)]
        ranges_cd2 = [(r.l2[i], r.r2[i]) for i in range(r.nblock)]
        cand_map[key] = {
            'score': r.score,
            'zscore': r.rmsd,  # parser maps rmsd field to zscore
            'rmsd': r.zmax,    # parser maps zmax field to rmsd
            'nblock': r.nblock,
            'ranges_cd1': ranges_cd1,
            'ranges_cd2': ranges_cd2,
        }

    all_keys = sorted(set(ref_map.keys()) | set(cand_map.keys()))
    matched = 0
    mismatched = 0
    missing = 0
    extra = 0

    score_tol = 1.0
    z_tol = 0.1

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
            if abs(ref.score - cand['score']) > score_tol:
                issues.append(f"score: ref={ref.score:.1f} cand={cand['score']:.1f}")
            if abs(ref.zscore - cand['zscore']) > z_tol:
                issues.append(f"zscore(rmsd): ref={ref.zscore:.1f} cand={cand['zscore']:.1f}")
            if ref.nblock != cand['nblock']:
                issues.append(f"nblock: ref={ref.nblock} cand={cand['nblock']}")
            if ref.ranges_cd1 != cand['ranges_cd1']:
                issues.append(f"cd1 ranges differ")
            if ref.ranges_cd2 != cand['ranges_cd2']:
                issues.append(f"cd2 ranges differ")

            if issues:
                print(f"  MISMATCH: {pair_label}: {'; '.join(issues)}")
                mismatched += 1
            else:
                print(f"  MATCH:   {pair_label} (score={cand['score']:.1f}, nblock={cand['nblock']})")
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
