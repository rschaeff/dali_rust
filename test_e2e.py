#!/usr/bin/env python3
"""
End-to-end pipeline validation against reference DaliLite output.

For each reference pair, runs the full pipeline (both wolf and parsi paths,
both directions) and compares DCCP output against reference.
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from pipeline import compare_pair
from validation.parsers import parse_dccp

DAT_DIR = Path('validation/fixtures/structures')
E2E_DIR = Path('validation/fixtures/e2e')

PAIRS = [
    ('101mA', '1a00A'),
    ('101mA', '1binA'),
    ('1a87A', '1allA'),
]

# parse_dccp stores zmax in the .rmsd field and rmsd in the .zscore field
# due to positional parsing. Our DCCPResult has the fields named correctly.
SCORE_TOL = 2.0


def compare_dccp_blocks(cand_results, ref_blocks, cd1, cd2, label):
    """Compare candidate DCCPResults against reference DCCPBlocks.

    Groups by direction (cd1->cd2 vs cd2->cd1) and compares scores.

    Returns:
        (matched, total, details)
    """
    details = []
    matched = 0
    total = 0

    for direction_cd1, direction_cd2, dir_label in [
        (cd1, cd2, f'{cd1}->{cd2}'),
        (cd2, cd1, f'{cd2}->{cd1}'),
    ]:
        ref_dir = [b for b in ref_blocks
                   if b.cd1 == direction_cd1 and b.cd2 == direction_cd2]
        cand_dir = [r for r in cand_results
                    if r.cd1 == direction_cd1 and r.cd2 == direction_cd2]

        if not ref_dir and not cand_dir:
            continue

        total += max(len(ref_dir), len(cand_dir))

        if not ref_dir:
            details.append(f"    {dir_label}: {len(cand_dir)} EXTRA results (ref has 0)")
            continue
        if not cand_dir:
            details.append(f"    {dir_label}: MISSING (ref has {len(ref_dir)})")
            continue

        # Compare count
        if len(cand_dir) != len(ref_dir):
            details.append(f"    {dir_label}: count mismatch cand={len(cand_dir)} "
                           f"ref={len(ref_dir)}")

        # Compare scores: ref.rmsd is actually zmax due to parser field order
        for i in range(min(len(cand_dir), len(ref_dir))):
            r = cand_dir[i]
            ref = ref_dir[i]
            score_diff = abs(r.score - ref.score)
            # ref.rmsd holds zmax, ref.zscore holds rmsd
            z_diff = abs(r.zmax - ref.rmsd)

            ok_score = score_diff <= SCORE_TOL
            ok_z = z_diff <= 0.5

            if ok_score:
                matched += 1
                details.append(
                    f"    {dir_label}[{i}]: MATCH score={r.score:.1f} "
                    f"(ref={ref.score:.1f} diff={score_diff:.1f}) "
                    f"z={r.zmax:.1f} (ref={ref.rmsd:.1f} diff={z_diff:.1f})")
            else:
                details.append(
                    f"    {dir_label}[{i}]: DIFF  score={r.score:.1f} "
                    f"(ref={ref.score:.1f} diff={score_diff:.1f}) "
                    f"z={r.zmax:.1f} (ref={ref.rmsd:.1f} diff={z_diff:.1f})")

    return matched, total, details


def test_pair(cd1, cd2):
    """Run full pipeline on one pair and compare against reference."""
    pair_dir = E2E_DIR / f'{cd1}_vs_{cd2}'
    ref_path = pair_dir / f'{cd1}.dccp'

    ref_blocks = parse_dccp(str(ref_path))
    # Filter self-comparisons from reference (shouldn't be any, but just in case)
    ref_blocks = [b for b in ref_blocks if b.cd1 != b.cd2]

    print(f"\n  --- {cd1} vs {cd2} ---")
    print(f"  Reference: {len(ref_blocks)} DCCP blocks")

    t0 = time.time()
    results = compare_pair(cd1, cd2, str(DAT_DIR))
    elapsed = time.time() - t0
    print(f"  Pipeline:  {len(results)} DCCP results ({elapsed:.0f}s)")

    matched, total, details = compare_dccp_blocks(
        results, ref_blocks, cd1, cd2, f'{cd1}_vs_{cd2}')

    for d in details:
        print(d)

    return matched, total


def main():
    print("=" * 60)
    print("End-to-End Pipeline Validation")
    print("=" * 60)

    total_matched = 0
    total_count = 0

    for cd1, cd2 in PAIRS:
        matched, total = test_pair(cd1, cd2)
        total_matched += matched
        total_count += total

    print(f"\n{'=' * 60}")
    print(f"TOTAL: {total_matched}/{total_count} direction-entries matched "
          f"(score tolerance={SCORE_TOL})")
    print(f"{'=' * 60}")

    return total_matched == total_count


if __name__ == '__main__':
    ok = main()
    sys.exit(0 if ok else 1)
