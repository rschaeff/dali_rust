#!/usr/bin/env python3
"""
Batch comparison: Rust iterative_search vs Fortran iterative DALI.

Compares domain detection at the biological level:
  - How many domains (masking rounds) does each method find?
  - For each Fortran domain, does Rust detect a domain covering the same region?
  - For shared domains, do z-scores agree?

Fortran iterative DALI: tests ALL templates per round, masks best hit, repeats.
  Each round = one domain. Multiple templates per round = alternatives for that domain.

Rust iterative_search: finds globally best template per round, masks, repeats.
  Each round = one domain with one template assignment.

Accuracy = do they find the same domains (query regions)?

Usage:
    python3 test_batch_comparison.py manifest.json
"""

import json
import sys
import time
import logging
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, '/home/rschaeff/dev/dpam_c2')
logging.basicConfig(level=logging.WARNING)

ECOD70_PDB_DIR = Path('/home/rschaeff/data/dpam_reference/ecod_data/ECOD70')

# Parse args: manifest [--dat-dir DIR]
import argparse
_parser = argparse.ArgumentParser()
_parser.add_argument('manifest', nargs='?', default='/tmp/batch_test_manifest.json')
_parser.add_argument('--dat-dir', default=None, help='Shared .dat cache dir (reused across proteins)')
_args = _parser.parse_args()
MANIFEST = _args.manifest
DAT_DIR = _args.dat_dir


def parse_fortran_hits(hits_file):
    """Parse Fortran _iterativdDali_hits file.

    Returns list of dicts with template, round, zscore, query_residues (set).
    """
    hits = []
    current = None
    with open(hits_file) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if current is not None:
                    hits.append(current)
                parts = line[1:].split('\t')
                # >TEMPLATE_ROUND\tZSCORE\tNRES1\tNRES2\tIDENTITY
                template_round = parts[0]
                template = template_round.rsplit('_', 1)[0]
                round_num = int(template_round.rsplit('_', 1)[1])
                current = {
                    'template': template,
                    'round': round_num,
                    'zscore': float(parts[1]),
                    'nres_query': int(parts[2]),
                    'nres_template': int(parts[3]),
                    'query_residues': set(),
                }
            elif current is not None:
                if line.startswith('rotation') or line.startswith('translation'):
                    continue
                parts = line.split('\t')
                if len(parts) == 2:
                    q_idx = int(parts[0])
                    current['query_residues'].add(q_idx)
        if current is not None:
            hits.append(current)
    return hits


def fortran_domains(hits):
    """Group Fortran hits into domains (one per round).

    Returns dict of round → {best_template, best_zscore, all_residues (union),
                              n_templates, z_range}
    """
    rounds = defaultdict(list)
    for h in hits:
        rounds[h['round']].append(h)

    domains = {}
    for r, round_hits in sorted(rounds.items()):
        # Best hit in this round (by z-score)
        best = max(round_hits, key=lambda x: x['zscore'])
        # Union of all query residues across templates in this round
        all_residues = set()
        for h in round_hits:
            all_residues.update(h['query_residues'])
        domains[r] = {
            'round': r,
            'best_template': best['template'],
            'best_zscore': best['zscore'],
            'query_residues': all_residues,
            'n_templates': len(round_hits),
            'templates': {h['template']: h['zscore'] for h in round_hits},
        }
    return domains


def rust_domains(results):
    """Convert Rust iterative_search results to domain list.

    Each result = one domain (one round of iterative search).
    Returns dict of round → {template, zscore, query_residues}
    """
    domains = {}
    for i, r in enumerate(results):
        template, zscore, n_aligned, qlen, alignments = r[0], r[1], r[2], r[3], r[4]
        query_residues = set(a[0] for a in alignments)
        domains[i + 1] = {
            'round': i + 1,
            'template': template,
            'zscore': zscore,
            'n_aligned': n_aligned,
            'query_residues': query_residues,
        }
    return domains


def region_overlap(set_a, set_b):
    """Jaccard-like overlap: intersection / min(|A|, |B|)."""
    if not set_a or not set_b:
        return 0.0
    return len(set_a & set_b) / min(len(set_a), len(set_b))


def match_domains(fortran_doms, rust_doms, min_overlap=0.5):
    """Match Fortran and Rust domains by query region overlap.

    Returns:
      matched: list of (fortran_round, rust_round, overlap, template_match)
      fortran_only: list of fortran_rounds with no Rust match
      rust_only: list of rust_rounds with no Fortran match
    """
    # Greedy matching: best overlap first
    pairs = []
    for fr, fdom in fortran_doms.items():
        for rr, rdom in rust_doms.items():
            ovl = region_overlap(fdom['query_residues'], rdom['query_residues'])
            if ovl >= min_overlap:
                pairs.append((fr, rr, ovl, rdom['template'] in fdom['templates']))

    # Sort by overlap descending, then greedily assign
    pairs.sort(key=lambda x: -x[2])
    used_f, used_r = set(), set()
    matched = []
    for fr, rr, ovl, tmatch in pairs:
        if fr not in used_f and rr not in used_r:
            matched.append((fr, rr, ovl, tmatch))
            used_f.add(fr)
            used_r.add(rr)

    fortran_only = [r for r in fortran_doms if r not in used_f]
    rust_only = [r for r in rust_doms if r not in used_r]
    return matched, fortran_only, rust_only


def run_batch_test():
    import dali as dali_rust
    from dpam.tools.dali import RustDALI

    with open(MANIFEST) as f:
        proteins = json.load(f)

    print(f"Batch accuracy test: {len(proteins)} proteins")
    print(f"Manifest: {MANIFEST}")
    print(f"{'='*90}")

    summary = {
        'total': len(proteins),
        'success': 0,
        'error': 0,
        'no_fortran': 0,
        'no_rust': 0,
    }

    # Aggregate counters
    total_fortran_domains = 0
    total_rust_domains = 0
    total_matched = 0
    total_fortran_only = 0
    total_rust_only = 0
    total_template_match = 0  # matched domains where Rust template ∈ Fortran templates
    all_zscore_pairs = []  # (rust_z, fortran_best_z) for matched domains
    all_fortran_missed = []  # (template, zscore, n_templates) for unmatched Fortran domains
    timings = []
    results_detail = []

    for i, protein in enumerate(proteins):
        name = protein['name']
        short_name = name[:50]
        pdb_file = Path(protein['pdb'])
        hits4dali = Path(protein['hits4dali'])
        fortran_file = Path(protein['fortran_hits'])

        # Read candidate templates
        edomains = [line.strip() for line in open(hits4dali) if line.strip()]

        # Parse Fortran reference into domains
        fortran_hits = parse_fortran_hits(fortran_file)
        f_doms = fortran_domains(fortran_hits)
        n_f = len(f_doms)

        if n_f == 0:
            summary['no_fortran'] += 1

        # Run Rust iterative search
        dali = RustDALI(dat_dir=DAT_DIR)
        try:
            t0 = time.time()
            rust_results = dali.batch_search(
                pdb_file, edomains,
                pdb_dir=ECOD70_PDB_DIR,
                min_aligned=20, min_zscore=2.0,
                gap_tolerance=5, max_rounds=len(edomains),
            )
            elapsed = time.time() - t0
            timings.append(elapsed)
        except Exception as e:
            print(f"[{i+1:3d}/{len(proteins)}] {short_name:50s}  ERROR: {e}")
            summary['error'] += 1
            results_detail.append({
                'name': name, 'error': str(e),
                'n_fortran_domains': n_f, 'n_rust_domains': 0,
            })
            continue

        summary['success'] += 1
        r_doms = rust_domains(rust_results)
        n_r = len(r_doms)

        if n_r == 0:
            summary['no_rust'] += 1

        # Match domains
        matched, f_only, r_only = match_domains(f_doms, r_doms)

        total_fortran_domains += n_f
        total_rust_domains += n_r
        total_matched += len(matched)
        total_fortran_only += len(f_only)
        total_rust_only += len(r_only)

        for fr, rr, ovl, tmatch in matched:
            if tmatch:
                total_template_match += 1
            all_zscore_pairs.append((r_doms[rr]['zscore'], f_doms[fr]['best_zscore']))

        for fr in f_only:
            fd = f_doms[fr]
            all_fortran_missed.append((fd['best_template'], fd['best_zscore'], fd['n_templates']))

        # Per-protein output
        status = "OK" if len(f_only) == 0 and len(r_only) == 0 else "DIFF"
        print(f"[{i+1:3d}/{len(proteins)}] {short_name:50s}  "
              f"F:{n_f:2d} R:{n_r:2d}  "
              f"match:{len(matched):2d} Fonly:{len(f_only):2d} Ronly:{len(r_only):2d}  "
              f"{elapsed:5.1f}s  {status}")

        for fr, rr, ovl, tmatch in matched:
            fd, rd = f_doms[fr], r_doms[rr]
            tag = "TMATCH" if tmatch else "REGION"
            print(f"    {tag:7s}  F{fr}({fd['best_template']}  z={fd['best_zscore']:5.1f}  n={fd['n_templates']:3d}) "
                  f"↔ R{rr}({rd['template']}  z={rd['zscore']:5.1f})  ovl={ovl:.2f}")

        for fr in f_only:
            fd = f_doms[fr]
            print(f"    MISSED   F{fr}({fd['best_template']}  z={fd['best_zscore']:5.1f}  n={fd['n_templates']:3d})")

        for rr in r_only:
            rd = r_doms[rr]
            print(f"    NOVEL    R{rr}({rd['template']}  z={rd['zscore']:5.1f})")

        results_detail.append({
            'name': name,
            'n_templates': len(edomains),
            'n_fortran_domains': n_f,
            'n_rust_domains': n_r,
            'n_matched': len(matched),
            'n_fortran_only': len(f_only),
            'n_rust_only': len(r_only),
            'elapsed': elapsed,
            'matched': [(fr, rr, ovl, tmatch) for fr, rr, ovl, tmatch in matched],
            'fortran_only': [(fr, f_doms[fr]['best_template'], f_doms[fr]['best_zscore'])
                            for fr in f_only],
            'rust_only': [(rr, r_doms[rr]['template'], r_doms[rr]['zscore'])
                         for rr in r_only],
        })

    # ==================== SUMMARY ====================
    print(f"\n{'='*90}")
    print(f"BATCH ACCURACY SUMMARY")
    print(f"{'='*90}")
    print(f"Proteins tested:    {summary['total']}")
    print(f"Successful:         {summary['success']}")
    print(f"Errors:             {summary['error']}")
    print(f"No Fortran domains: {summary['no_fortran']}")
    print(f"No Rust domains:    {summary['no_rust']}")

    print(f"\nDomain detection:")
    print(f"  Fortran domains:    {total_fortran_domains}")
    print(f"  Rust domains:       {total_rust_domains}")
    print(f"  Matched:            {total_matched} ({total_matched/max(total_fortran_domains,1)*100:.1f}% recall)")
    print(f"  Fortran-only:       {total_fortran_only} (missed by Rust)")
    print(f"  Rust-only:          {total_rust_only} (novel Rust detections)")
    if total_matched > 0:
        print(f"  Template agreement: {total_template_match}/{total_matched} "
              f"({total_template_match/total_matched*100:.1f}% of matched domains)")

    # Recall stratified by Fortran best z-score
    z_thresholds = [2.0, 5.0, 10.0, 15.0, 20.0]
    print(f"\nRecall by Fortran best z-score:")
    print(f"  {'z≥':>6s}  {'F_doms':>6s}  {'matched':>7s}  {'recall':>7s}")
    # Collect all fortran domain z-scores
    all_f_zscores = []  # (zscore, was_matched)
    for r in results_detail:
        if 'error' in r:
            continue
        for fr, rr, ovl, tmatch in r['matched']:
            # Find fortran z for this round
            all_f_zscores.append((r['matched'][r['matched'].index((fr, rr, ovl, tmatch))], True))
        for fr, ft, fz in r['fortran_only']:
            all_f_zscores.append((fz, False))

    # Simpler: just collect from matched + missed with z-scores
    matched_fz = [fz for rz, fz in all_zscore_pairs]
    missed_fz = [fz for _, fz, _ in all_fortran_missed]
    all_fz = [(z, True) for z in matched_fz] + [(z, False) for z in missed_fz]

    for z_thresh in z_thresholds:
        above = [(z, m) for z, m in all_fz if z >= z_thresh]
        found = sum(1 for _, m in above if m)
        pct = found / max(len(above), 1) * 100
        print(f"  {z_thresh:6.1f}  {len(above):6d}  {found:7d}  {pct:6.1f}%")

    # Z-score agreement
    if all_zscore_pairs:
        diffs = [rz - fz for rz, fz in all_zscore_pairs]
        diffs.sort()
        n = len(diffs)
        mean_d = sum(diffs) / n
        print(f"\nZ-score agreement (matched domains, n={n}):")
        print(f"  Mean diff (Rust-Fortran): {mean_d:+.2f}")
        print(f"  Median diff:              {diffs[n//2]:+.2f}")
        print(f"  Std dev:                  {(sum((d - mean_d)**2 for d in diffs)/n)**0.5:.2f}")
        print(f"  Within ±1:  {sum(1 for d in diffs if abs(d) <= 1)/n*100:.1f}%")
        print(f"  Within ±2:  {sum(1 for d in diffs if abs(d) <= 2)/n*100:.1f}%")
        print(f"  Within ±5:  {sum(1 for d in diffs if abs(d) <= 5)/n*100:.1f}%")

    # Missed domain analysis
    if all_fortran_missed:
        missed_z = sorted([fz for _, fz, _ in all_fortran_missed], reverse=True)
        print(f"\nMissed Fortran domains (n={len(missed_z)}):")
        print(f"  Z-score range: {missed_z[-1]:.1f} – {missed_z[0]:.1f}")
        print(f"  Z ≥ 10: {sum(1 for z in missed_z if z >= 10)}")
        print(f"  Z 5-10: {sum(1 for z in missed_z if 5 <= z < 10)}")
        print(f"  Z 2-5:  {sum(1 for z in missed_z if 2 <= z < 5)}")

    # Timing
    if timings:
        timings.sort()
        n = len(timings)
        print(f"\nTiming:")
        print(f"  Total:  {sum(timings):.1f}s ({sum(timings)/3600:.1f}h)")
        print(f"  Mean:   {sum(timings)/n:.1f}s per protein")
        print(f"  Median: {timings[n//2]:.1f}s")
        print(f"  P95:    {timings[int(n*0.95)]:.1f}s")
        print(f"  Max:    {timings[-1]:.1f}s")

    # Save detailed results
    output_path = MANIFEST.replace('.json', '_results.json')
    with open(output_path, 'w') as f:
        json.dump(results_detail, f, indent=2)
    print(f"\nDetailed results: {output_path}")


if __name__ == '__main__':
    run_batch_test()
