#!/usr/bin/env python3
"""Test PARSI and FILTER95 modules against 18-structure expanded corpus ground truth."""

import sys
import time
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).parent))

STRUCTURES = [
    '1a7sA', '1aiwA', '1a25A', '1a12A', '1f3uA', '1a04A',
    '1bbhA', '1b3qA', '1a17A', '1miwA', '1aopA', '1a8lA',
    '1a6qA', '1a06A', '1b3oB', '1a0cA', '1a4iA', '1bcoA',
]

GT_DIR = Path('validation/corpus_expansion/ground_truth_18')
DAT_DIR = Path('validation/corpus_expansion/dat')

# Also try DaliLite DAT dir for structures
STRUCTURE_DAT_DIR = None
for candidate in [DAT_DIR, Path('/home/rschaeff/src/Dali_v5/DaliLite.v5/DAT')]:
    if candidate.exists():
        STRUCTURE_DAT_DIR = candidate
        break

if STRUCTURE_DAT_DIR is None:
    print("ERROR: Could not find .dat files directory")
    sys.exit(1)

print(f"Using .dat files from: {STRUCTURE_DAT_DIR}")


def parse_refine_lines(lines):
    """Parse PARSI refine output lines."""
    parsed = []
    for line in lines:
        line = line.strip()
        if not line or 'refine' not in line:
            continue
        parts = line.split()
        cd1cd2 = parts[0].replace('refine', '')
        idom = int(parts[1])
        score = int(parts[2])
        nseg = int(parts[3])
        ranges = [int(x) for x in parts[4:]]
        parsed.append({
            'cd1cd2': cd1cd2, 'idom': idom, 'score': score,
            'nseg': nseg, 'ranges': ranges,
        })
    return parsed


def parse_filter95_output(lines):
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
        parsed[key] = {'zscore': zscore, 'score': score, 'nres': nres}
    return parsed


def test_filter95():
    """Run FILTER95 on 18-structure corpus and compare to ground truth."""
    from filter95.filter95 import run_filter95

    print("=" * 70)
    print("FILTER95 Module - Expanded corpus (18 structures)")
    print("=" * 70)

    gt_input = GT_DIR / 'filter95' / 'filter95_input.txt'
    gt_output = GT_DIR / 'filter95' / 'filter95_output.txt'

    with open(gt_input) as f:
        input_lines = f.readlines()
    print(f"  Input: {len(input_lines)} refine lines")

    with open(gt_output) as f:
        ref_text = f.read()
    ref_lines = [l.strip() for l in ref_text.strip().split('\n') if l.strip()]
    print(f"  Reference: {len(ref_lines)} output lines")

    t0 = time.time()
    result_lines = run_filter95(input_lines, str(STRUCTURE_DAT_DIR) + '/')
    elapsed = time.time() - t0
    print(f"  Python: {len(result_lines)} output lines ({elapsed:.1f}s)")

    ref_parsed = parse_filter95_output(ref_lines)
    cand_parsed = parse_filter95_output(result_lines)

    matched = 0
    mismatched = 0
    for key, ref in ref_parsed.items():
        cand = cand_parsed.get(key)
        if cand is None:
            mismatched += 1
            if mismatched <= 10:
                print(f"  MISSING: {key} zscore={ref['zscore']:.3f} score={ref['score']}")
            continue

        z_ok = abs(ref['zscore'] - cand['zscore']) < 0.5
        s_ok = abs(ref['score'] - cand['score']) <= max(abs(ref['score']) * 0.01, 100)

        if z_ok and s_ok:
            matched += 1
        else:
            mismatched += 1
            if mismatched <= 10:
                print(f"  MISMATCH: {key}")
                print(f"    ref:  z={ref['zscore']:.3f} score={ref['score']}")
                print(f"    cand: z={cand['zscore']:.3f} score={cand['score']}")

    extra = len(cand_parsed) - len(ref_parsed)
    print(f"\n  FILTER95: {matched}/{len(ref_parsed)} matched, "
          f"{mismatched} missing/mismatched, {extra} extra")
    return matched, mismatched


def test_parsi():
    """Run PARSI on 18-structure corpus and compare to ground truth."""
    from parsi.parsi import dowork_parsi, init_parsi

    print("\n" + "=" * 70)
    print("PARSI Module - Expanded corpus (18 structures)")
    print("=" * 70)

    gt_output = GT_DIR / 'parsi' / 'parsi_output.txt'
    with open(gt_output) as f:
        ref_lines = [l.strip() for l in f if l.strip()]
    ref_parsed = parse_refine_lines(ref_lines)
    print(f"  Reference: {len(ref_parsed)} refine lines", flush=True)

    t0 = time.time()
    init_parsi()
    all_lines = []
    cd1_cache = {}
    dat_dir = str(STRUCTURE_DAT_DIR) + '/'
    total_pairs = len(STRUCTURES) * len(STRUCTURES)
    pair_num = 0
    for cd1 in STRUCTURES:
        for cd2 in STRUCTURES:
            pair_num += 1
            t1 = time.time()
            lines = dowork_parsi(cd1, cd2, dat_dir, lfirstonly=True,
                                 cd1_cache=cd1_cache)
            dt = time.time() - t1
            all_lines.extend(lines)
            elapsed = time.time() - t0
            print(f"  [{pair_num}/{total_pairs}] {cd1}-{cd2}: "
                  f"{len(lines)} lines ({dt:.1f}s, total {elapsed:.0f}s)",
                  flush=True)

    result_lines = all_lines
    elapsed = time.time() - t0
    cand_parsed = parse_refine_lines(result_lines)
    print(f"  Python: {len(cand_parsed)} refine lines ({elapsed:.1f}s)",
          flush=True)

    # Group by (cd1cd2, idom)
    ref_map = defaultdict(list)
    for r in ref_parsed:
        ref_map[(r['cd1cd2'], r['idom'])].append(r)

    cand_map = defaultdict(list)
    for c in cand_parsed:
        cand_map[(c['cd1cd2'], c['idom'])].append(c)

    matched = 0
    mismatched = 0

    for key in sorted(ref_map.keys()):
        ref_group = ref_map[key]
        cand_group = cand_map.get(key, [])

        for ref in ref_group:
            found = False
            # Try exact match first
            for cand in cand_group:
                if (ref['score'] == cand['score'] and
                        ref['nseg'] == cand['nseg'] and
                        ref['ranges'] == cand['ranges']):
                    found = True
                    break
            if not found:
                # Try with tolerance
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

    total = matched + mismatched
    pct = 100.0 * matched / total if total else 0
    print(f"\n  PARSI: {matched}/{total} ({pct:.1f}%) matched, {mismatched} missing")
    return matched, mismatched


if __name__ == '__main__':
    if '--filter95' in sys.argv:
        matched, mismatched = test_filter95()
        sys.exit(0 if mismatched == 0 else 1)
    elif '--parsi' in sys.argv:
        matched, mismatched = test_parsi()
        sys.exit(0 if mismatched == 0 else 1)
    elif '--both' in sys.argv:
        f_m, f_mm = test_filter95()
        p_m, p_mm = test_parsi()
        sys.exit(0 if f_mm == 0 and p_mm == 0 else 1)
    else:
        # Default: filter95 only (PARSI takes hours)
        test_filter95()
