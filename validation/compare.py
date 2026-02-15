#!/usr/bin/env python3
"""
Compare reimplementation output against DaliLite ground truth fixtures.

Provides both stage-level and end-to-end comparison with configurable
tolerances for floating-point scores.

Usage:
    # Compare a new WOLF implementation's output
    python compare.py wolf --reference fixtures/wolf/wolf_output.txt \
                           --candidate /path/to/new_wolf_output.txt

    # Compare DP Z-scores
    python compare.py dp --reference fixtures/dp/wolf/dp_output.txt \
                         --candidate /path/to/new_dp_output.txt

    # Compare end-to-end DCCP results
    python compare.py dccp --reference fixtures/e2e/2nrmA_vs_101mA/2nrmA.dccp \
                           --candidate /path/to/new_2nrmA.dccp

    # Run full validation suite
    python compare.py suite --fixtures-dir fixtures/ --candidate-dir /path/to/new/output/
"""

import argparse
import json
import sys
from pathlib import Path

from parsers import (
    parse_wolfitz,
    parse_dccp,
    parse_refine_lines,
    parse_filter95,
)


class ComparisonResult:
    def __init__(self, stage, label=''):
        self.stage = stage
        self.label = label
        self.total = 0
        self.matched = 0
        self.mismatched = 0
        self.missing_in_candidate = 0
        self.extra_in_candidate = 0
        self.details = []

    def add_match(self, key, detail=''):
        self.total += 1
        self.matched += 1
        if detail:
            self.details.append(('MATCH', key, detail))

    def add_mismatch(self, key, detail):
        self.total += 1
        self.mismatched += 1
        self.details.append(('MISMATCH', key, detail))

    def add_missing(self, key, detail=''):
        self.total += 1
        self.missing_in_candidate += 1
        self.details.append(('MISSING', key, detail))

    def add_extra(self, key, detail=''):
        self.total += 1
        self.extra_in_candidate += 1
        self.details.append(('EXTRA', key, detail))

    @property
    def passed(self):
        return self.mismatched == 0 and self.missing_in_candidate == 0

    def summary(self):
        status = 'PASS' if self.passed else 'FAIL'
        label = f' ({self.label})' if self.label else ''
        lines = [f'[{status}] {self.stage}{label}: '
                 f'{self.matched}/{self.total} matched, '
                 f'{self.mismatched} mismatched, '
                 f'{self.missing_in_candidate} missing, '
                 f'{self.extra_in_candidate} extra']
        for dtype, key, detail in self.details:
            if dtype != 'MATCH':
                lines.append(f'  {dtype}: {key} - {detail}')
        return '\n'.join(lines)


def compare_wolf(ref_path, cand_path, **kwargs):
    """Compare WOLF outputs: same pairs found, same nblock counts."""
    result = ComparisonResult('WOLF')

    ref_hits = parse_wolfitz(ref_path)
    cand_hits = parse_wolfitz(cand_path)

    # Index by (cd1, cd2)
    ref_map = {}
    for h in ref_hits:
        key = (h.cd1, h.cd2)
        ref_map[key] = h

    cand_map = {}
    for h in cand_hits:
        key = (h.cd1, h.cd2)
        cand_map[key] = h

    all_keys = set(ref_map.keys()) | set(cand_map.keys())
    for key in sorted(all_keys):
        cd1, cd2 = key
        pair_label = f'{cd1}-{cd2}'
        if key not in cand_map:
            result.add_missing(pair_label, 'pair not found in candidate')
        elif key not in ref_map:
            result.add_extra(pair_label, 'pair not in reference')
        else:
            ref_h = ref_map[key]
            cand_h = cand_map[key]
            if ref_h.nblock != cand_h.nblock:
                result.add_mismatch(pair_label,
                                    f'nblock: ref={ref_h.nblock} cand={cand_h.nblock}')
            elif ref_h.raw_values != cand_h.raw_values:
                result.add_mismatch(pair_label, 'alignment values differ')
            else:
                result.add_match(pair_label)

    return result


def compare_dccp(ref_path, cand_path, score_tol=1.0, z_tol=0.1, **kwargs):
    """Compare DCCP outputs: scores within tolerance, same alignment blocks."""
    result = ComparisonResult('DCCP')

    ref_blocks = parse_dccp(ref_path)
    cand_blocks = parse_dccp(cand_path)

    # Group by (cd1, cd2) — may have multiple blocks per pair
    def group_by_pair(blocks):
        groups = {}
        for b in blocks:
            key = (b.cd1, b.cd2)
            groups.setdefault(key, []).append(b)
        return groups

    ref_groups = group_by_pair(ref_blocks)
    cand_groups = group_by_pair(cand_blocks)

    all_keys = set(ref_groups.keys()) | set(cand_groups.keys())
    for key in sorted(all_keys):
        cd1, cd2 = key
        pair_label = f'{cd1}-{cd2}'

        if key not in cand_groups:
            result.add_missing(pair_label,
                               f'{len(ref_groups[key])} blocks in reference')
            continue
        if key not in ref_groups:
            result.add_extra(pair_label,
                             f'{len(cand_groups[key])} blocks in candidate')
            continue

        # Compare best block per pair (highest score)
        ref_best = max(ref_groups[key], key=lambda b: b.score)
        cand_best = max(cand_groups[key], key=lambda b: b.score)

        issues = []
        if abs(ref_best.score - cand_best.score) > score_tol:
            issues.append(
                f'score: ref={ref_best.score:.1f} cand={cand_best.score:.1f}')
        if abs(ref_best.zscore - cand_best.zscore) > z_tol:
            issues.append(
                f'zscore: ref={ref_best.zscore:.1f} cand={cand_best.zscore:.1f}')
        if ref_best.nblock != cand_best.nblock:
            issues.append(
                f'nblock: ref={ref_best.nblock} cand={cand_best.nblock}')
        if ref_best.ranges_cd1 != cand_best.ranges_cd1:
            issues.append('cd1 ranges differ')
        if ref_best.ranges_cd2 != cand_best.ranges_cd2:
            issues.append('cd2 ranges differ')

        if issues:
            result.add_mismatch(pair_label, '; '.join(issues))
        else:
            result.add_match(pair_label)

    return result


def compare_refine(ref_path, cand_path, score_tol=100, **kwargs):
    """Compare refine line outputs from PARSI/FILTER95."""
    result = ComparisonResult('REFINE')

    ref_lines = parse_refine_lines(ref_path)
    cand_lines = parse_refine_lines(cand_path)

    # Group by (cd1cd2, idom)
    def group(lines):
        groups = {}
        for r in lines:
            key = (r.cd1cd2, r.idom)
            groups.setdefault(key, []).append(r)
        return groups

    ref_groups = group(ref_lines)
    cand_groups = group(cand_lines)

    all_keys = set(ref_groups.keys()) | set(cand_groups.keys())
    for key in sorted(all_keys):
        label = f'{key[0]}-dom{key[1]}'
        if key not in cand_groups:
            result.add_missing(label)
        elif key not in ref_groups:
            result.add_extra(label)
        else:
            # Compare best score per group
            ref_best = max(ref_groups[key], key=lambda r: r.score)
            cand_best = max(cand_groups[key], key=lambda r: r.score)

            if abs(ref_best.score - cand_best.score) > score_tol:
                result.add_mismatch(label,
                    f'score: ref={ref_best.score} cand={cand_best.score}')
            elif ref_best.ranges_cd1 != cand_best.ranges_cd1:
                result.add_mismatch(label, 'cd1 ranges differ')
            elif ref_best.ranges_cd2 != cand_best.ranges_cd2:
                result.add_mismatch(label, 'cd2 ranges differ')
            else:
                result.add_match(label)

    return result


def compare_filter95(ref_path, cand_path, z_tol=0.5, score_tol=100000, **kwargs):
    """Compare FILTER95 outputs: same pairs, Z-scores within tolerance."""
    result = ComparisonResult('FILTER95')

    ref_lines = parse_filter95(ref_path)
    cand_lines = parse_filter95(cand_path)

    # Group by (cd1, cd2) - take best Z-score per pair
    def best_per_pair(lines):
        best = {}
        for line in lines:
            key = (line.cd1, line.cd2)
            if key not in best or line.zscore > best[key].zscore:
                best[key] = line
        return best

    ref_best = best_per_pair(ref_lines)
    cand_best = best_per_pair(cand_lines)

    all_keys = set(ref_best.keys()) | set(cand_best.keys())
    for key in sorted(all_keys):
        cd1, cd2 = key
        pair_label = f'{cd1}-{cd2}'
        if key not in cand_best:
            result.add_missing(pair_label)
        elif key not in ref_best:
            result.add_extra(pair_label)
        else:
            ref_l = ref_best[key]
            cand_l = cand_best[key]
            issues = []
            if abs(ref_l.zscore - cand_l.zscore) > z_tol:
                issues.append(
                    f'zscore: ref={ref_l.zscore:.1f} cand={cand_l.zscore:.1f}')
            if abs(ref_l.score - cand_l.score) > score_tol:
                issues.append(
                    f'score: ref={ref_l.score} cand={cand_l.score}')
            if ref_l.nres != cand_l.nres:
                issues.append(
                    f'nres: ref={ref_l.nres} cand={cand_l.nres}')
            if issues:
                result.add_mismatch(pair_label, '; '.join(issues))
            else:
                result.add_match(pair_label)

    return result


def compare_rankings(ref_path, cand_path, **kwargs):
    """Compare hit rankings: same ordering of targets by Z-score."""
    result = ComparisonResult('RANKING')

    ref_blocks = parse_dccp(ref_path)
    cand_blocks = parse_dccp(cand_path)

    def best_per_pair(blocks):
        best = {}
        for b in blocks:
            key = (b.cd1, b.cd2)
            if key not in best or b.score > best[key].score:
                best[key] = b
        return best

    ref_best = best_per_pair(ref_blocks)
    cand_best = best_per_pair(cand_blocks)

    # Get unique cd1 values
    cd1_values = set(b.cd1 for b in ref_blocks) | set(b.cd1 for b in cand_blocks)

    for cd1 in sorted(cd1_values):
        ref_ranking = sorted(
            [(k, v) for k, v in ref_best.items() if k[0] == cd1],
            key=lambda x: -x[1].score)
        cand_ranking = sorted(
            [(k, v) for k, v in cand_best.items() if k[0] == cd1],
            key=lambda x: -x[1].score)

        ref_order = [k[1] for k, _ in ref_ranking]
        cand_order = [k[1] for k, _ in cand_ranking]

        if ref_order == cand_order:
            result.add_match(cd1, f'{len(ref_order)} targets')
        else:
            # Find first rank where they differ
            for i, (r, c) in enumerate(zip(ref_order, cand_order)):
                if r != c:
                    result.add_mismatch(cd1,
                        f'rank {i+1}: ref={r} cand={c}')
                    break

    return result


def run_suite(fixtures_dir, candidate_dir):
    """Run full validation suite against all captured fixtures."""
    results = []
    fixtures = Path(fixtures_dir)
    candidate = Path(candidate_dir)

    # WOLF comparison
    ref_wolf = fixtures / 'wolf' / 'wolf_output.txt'
    cand_wolf = candidate / 'wolf' / 'wolf_output.txt'
    if ref_wolf.exists() and cand_wolf.exists():
        results.append(compare_wolf(ref_wolf, cand_wolf))

    # DP comparisons
    for label in ('wolf', 'dalicon'):
        ref_dp = fixtures / 'dp' / label / 'dp_output.txt'
        cand_dp = candidate / 'dp' / label / 'dp_output.txt'
        if ref_dp.exists() and cand_dp.exists():
            r = compare_dccp(ref_dp, cand_dp)
            r.label = f'DP-{label}'
            results.append(r)

    # PARSI comparison
    ref_parsi = fixtures / 'parsi' / 'parsi_output.txt'
    cand_parsi = candidate / 'parsi' / 'parsi_output.txt'
    if ref_parsi.exists() and cand_parsi.exists():
        results.append(compare_refine(ref_parsi, cand_parsi))

    # FILTER95 comparison
    ref_f95 = fixtures / 'filter95' / 'filter95_output.txt'
    cand_f95 = candidate / 'filter95' / 'filter95_output.txt'
    if ref_f95.exists() and cand_f95.exists():
        results.append(compare_filter95(ref_f95, cand_f95))

    # E2E comparisons
    e2e_ref = fixtures / 'e2e'
    e2e_cand = candidate / 'e2e'
    if e2e_ref.exists() and e2e_cand.exists():
        for ref_case in sorted(e2e_ref.iterdir()):
            cand_case = e2e_cand / ref_case.name
            ref_dccp = list(ref_case.glob('*.dccp'))
            cand_dccp = list(cand_case.glob('*.dccp')) if cand_case.exists() else []
            if ref_dccp and cand_dccp:
                r = compare_dccp(ref_dccp[0], cand_dccp[0])
                r.label = ref_case.name
                results.append(r)
                # Also check rankings
                r2 = compare_rankings(ref_dccp[0], cand_dccp[0])
                r2.label = ref_case.name
                results.append(r2)

    return results


def main():
    parser = argparse.ArgumentParser(description='Compare against DaliLite ground truth')
    subparsers = parser.add_subparsers(dest='command')

    # Single-stage comparisons
    for stage in ('wolf', 'dp', 'dccp', 'refine', 'filter95', 'ranking'):
        sp = subparsers.add_parser(stage)
        sp.add_argument('--reference', required=True)
        sp.add_argument('--candidate', required=True)
        sp.add_argument('--score-tol', type=float, default=1.0)
        sp.add_argument('--z-tol', type=float, default=0.1)

    # Full suite
    sp = subparsers.add_parser('suite')
    sp.add_argument('--fixtures-dir', required=True)
    sp.add_argument('--candidate-dir', required=True)

    args = parser.parse_args()

    if args.command == 'suite':
        results = run_suite(args.fixtures_dir, args.candidate_dir)
        print(f'\n{"="*60}')
        print(f'VALIDATION SUITE RESULTS')
        print(f'{"="*60}\n')
        n_pass = sum(1 for r in results if r.passed)
        n_fail = sum(1 for r in results if not r.passed)
        for r in results:
            print(r.summary())
            print()
        print(f'{"="*60}')
        print(f'Total: {n_pass} passed, {n_fail} failed')
        sys.exit(0 if n_fail == 0 else 1)

    comparators = {
        'wolf': compare_wolf,
        'dp': compare_dccp,
        'dccp': compare_dccp,
        'refine': compare_refine,
        'filter95': compare_filter95,
        'ranking': compare_rankings,
    }

    func = comparators[args.command]
    result = func(args.reference, args.candidate,
                  score_tol=getattr(args, 'score_tol', 1.0),
                  z_tol=getattr(args, 'z_tol', 0.1))
    print(result.summary())
    sys.exit(0 if result.passed else 1)


if __name__ == '__main__':
    main()
