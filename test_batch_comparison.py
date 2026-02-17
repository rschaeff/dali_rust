#!/usr/bin/env python3
"""
Batch comparison: Rust batch_search vs Fortran iterative DALI.

Compares domain detection at the biological level:
  - How many domains does each method find?
  - Do they find the same ECOD templates?
  - Do the aligned regions overlap?

Usage:
    LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libgomp.so.1 python3 test_batch_comparison.py
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
MANIFEST = '/tmp/batch_test_manifest.json'


def parse_fortran_hits(hits_file):
    """Parse Fortran _iterativdDali_hits file.

    Returns list of dicts with template, zscore, n_aligned, query_residues (set).
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
                template = parts[0].rsplit('_', 1)[0]  # strip _N round suffix
                current = {
                    'template': template,
                    'zscore': float(parts[1]),
                    'nres_query': int(parts[2]),
                    'nres_template': int(parts[3]),
                    'query_residues': set(),
                    'alignments': [],
                }
            elif current is not None:
                if line.startswith('rotation') or line.startswith('translation'):
                    continue
                parts = line.split('\t')
                if len(parts) == 2:
                    q_idx, t_idx = int(parts[0]), int(parts[1])
                    current['query_residues'].add(q_idx)
                    current['alignments'].append((q_idx, t_idx))
        if current is not None:
            hits.append(current)
    return hits


def fortran_domains(hits):
    """Extract domain-level summary from Fortran hits.

    Groups per-template hits by covered query region.
    Returns list of (template, zscore, query_range_set) for the best hit
    per distinct template.
    """
    # Best hit per template (highest z-score)
    best = {}
    for h in hits:
        t = h['template']
        if t not in best or h['zscore'] > best[t]['zscore']:
            best[t] = h
    return best


def rust_domains(results):
    """Extract domain-level summary from Rust batch_search results.

    Returns dict of template → {zscore, query_residues, ...}
    """
    domains = {}
    for r in results:
        template, zscore, n_aligned, qlen, alignments = r[0], r[1], r[2], r[3], r[4]
        query_residues = set(a[0] for a in alignments)
        domains[template] = {
            'template': template,
            'zscore': zscore,
            'n_aligned': n_aligned,
            'query_residues': query_residues,
            'alignments': alignments,
        }
    return domains


def region_overlap(set_a, set_b):
    """Compute overlap fraction between two residue sets."""
    if not set_a or not set_b:
        return 0.0
    intersection = len(set_a & set_b)
    return intersection / min(len(set_a), len(set_b))


def find_matching_fortran_domain(rust_domain, fortran_best, min_overlap=0.5):
    """Find the Fortran hit that best matches a Rust domain by query coverage."""
    best_match = None
    best_overlap = 0.0
    for t, fhit in fortran_best.items():
        overlap = region_overlap(rust_domain['query_residues'], fhit['query_residues'])
        if overlap > best_overlap:
            best_overlap = overlap
            best_match = t
    if best_overlap >= min_overlap:
        return best_match, best_overlap
    return None, best_overlap


def run_batch_test():
    import dali as dali_rust
    from dpam.tools.dali import RustDALI

    with open(MANIFEST) as f:
        proteins = json.load(f)

    print(f"Batch test: {len(proteins)} proteins")
    print(f"{'='*80}")

    summary = {
        'total': len(proteins),
        'success': 0,
        'error': 0,
        'no_fortran_hits': 0,
        'no_rust_hits': 0,
        'domain_counts': [],      # (fortran_n, rust_n)
        'region_matches': [],     # fraction of Rust domains with Fortran overlap
        'template_matches': [],   # fraction of Rust domains with exact template match
        'timings': [],
    }

    results_detail = []

    for i, protein in enumerate(proteins):
        name = protein['name']
        short_name = name[:60]
        pdb_file = Path(protein['pdb'])
        hits4dali = Path(protein['hits4dali'])
        fortran_file = Path(protein['fortran_hits'])

        # Read inputs
        edomains = [line.strip() for line in open(hits4dali) if line.strip()]

        # Parse Fortran reference
        fortran_hits = parse_fortran_hits(fortran_file)
        fortran_best = fortran_domains(fortran_hits)
        n_fortran_domains = len(fortran_best)

        if n_fortran_domains == 0:
            summary['no_fortran_hits'] += 1

        # Run Rust batch search
        dali = RustDALI(dat_dir=None)
        try:
            t0 = time.time()
            rust_results = dali.batch_search(
                pdb_file, edomains,
                pdb_dir=ECOD70_PDB_DIR,
                min_aligned=20, min_zscore=2.0,
                gap_tolerance=5, max_rounds=len(edomains),
            )
            elapsed = time.time() - t0
            summary['timings'].append(elapsed)
        except Exception as e:
            print(f"[{i+1:3d}/{len(proteins)}] {short_name:60s}  ERROR: {e}")
            summary['error'] += 1
            continue

        summary['success'] += 1
        rust_doms = rust_domains(rust_results)
        n_rust_domains = len(rust_doms)

        if n_rust_domains == 0:
            summary['no_rust_hits'] += 1

        summary['domain_counts'].append((n_fortran_domains, n_rust_domains))

        # Compare: for each Rust domain, find best Fortran match
        n_region_match = 0
        n_template_match = 0
        match_details = []
        for template, rdom in rust_doms.items():
            # Exact template match
            if template in fortran_best:
                n_template_match += 1
                overlap = region_overlap(rdom['query_residues'],
                                        fortran_best[template]['query_residues'])
                match_details.append(('exact', template, rdom['zscore'],
                                     fortran_best[template]['zscore'], overlap))
            else:
                # Region-based match
                match_t, overlap = find_matching_fortran_domain(rdom, fortran_best)
                if match_t:
                    n_region_match += 1
                    match_details.append(('region', f"{template}~{match_t}",
                                         rdom['zscore'],
                                         fortran_best[match_t]['zscore'], overlap))
                else:
                    match_details.append(('novel', template, rdom['zscore'], 0, 0))

        if n_rust_domains > 0:
            summary['template_matches'].append(n_template_match / n_rust_domains)
            summary['region_matches'].append(
                (n_template_match + n_region_match) / n_rust_domains)
        else:
            summary['template_matches'].append(0)
            summary['region_matches'].append(0)

        # Print per-protein summary
        match_pct = ((n_template_match + n_region_match) / n_rust_domains * 100
                     if n_rust_domains > 0 else 0)
        print(f"[{i+1:3d}/{len(proteins)}] {short_name:60s}  "
              f"F:{n_fortran_domains:3d} R:{n_rust_domains:2d}  "
              f"match:{match_pct:5.1f}%  {elapsed:5.1f}s")

        for kind, desc, rz, fz, ovl in match_details:
            tag = {'exact': 'EXACT', 'region': 'REGION', 'novel': 'NOVEL'}[kind]
            print(f"         {tag:6s}  {desc:30s}  rz={rz:5.1f}  fz={fz:5.1f}  overlap={ovl:.2f}")

        results_detail.append({
            'name': name,
            'n_templates': len(edomains),
            'n_fortran_domains': n_fortran_domains,
            'n_rust_domains': n_rust_domains,
            'n_template_match': n_template_match,
            'n_region_match': n_region_match,
            'elapsed': elapsed,
            'matches': match_details,
        })

    # Print summary
    print(f"\n{'='*80}")
    print(f"BATCH TEST SUMMARY")
    print(f"{'='*80}")
    print(f"Proteins tested:    {summary['total']}")
    print(f"Successful:         {summary['success']}")
    print(f"Errors:             {summary['error']}")
    print(f"No Fortran hits:    {summary['no_fortran_hits']}")
    print(f"No Rust hits:       {summary['no_rust_hits']}")

    if summary['domain_counts']:
        fortran_counts = [d[0] for d in summary['domain_counts']]
        rust_counts = [d[1] for d in summary['domain_counts']]
        print(f"\nDomain counts:")
        print(f"  Fortran: mean={sum(fortran_counts)/len(fortran_counts):.1f}, "
              f"median={sorted(fortran_counts)[len(fortran_counts)//2]}")
        print(f"  Rust:    mean={sum(rust_counts)/len(rust_counts):.1f}, "
              f"median={sorted(rust_counts)[len(rust_counts)//2]}")

    if summary['region_matches']:
        region_matches = summary['region_matches']
        template_matches = summary['template_matches']
        print(f"\nDomain matching (per protein, Rust→Fortran):")
        print(f"  Exact template match:  mean={sum(template_matches)/len(template_matches)*100:.1f}%")
        print(f"  Region overlap match:  mean={sum(region_matches)/len(region_matches)*100:.1f}%")

    if summary['timings']:
        timings = summary['timings']
        print(f"\nTiming:")
        print(f"  Total:  {sum(timings):.1f}s")
        print(f"  Mean:   {sum(timings)/len(timings):.1f}s per protein")
        print(f"  Median: {sorted(timings)[len(timings)//2]:.1f}s")
        print(f"  Max:    {max(timings):.1f}s")

    # Save detailed results
    output_path = Path('/tmp/batch_test_results.json')
    with open(output_path, 'w') as f:
        # Convert sets/tuples for JSON serialization
        serializable = []
        for r in results_detail:
            r2 = dict(r)
            r2['matches'] = [(k, d, rz, fz, ovl) for k, d, rz, fz, ovl in r['matches']]
            serializable.append(r2)
        json.dump(serializable, f, indent=2)
    print(f"\nDetailed results: {output_path}")


if __name__ == '__main__':
    run_batch_test()
