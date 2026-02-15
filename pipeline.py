"""
End-to-end DaliLite pipeline orchestration.

Wires together:
  Path A: WOLF -> DP(zcut=2.0) -> DALICON(lfitz=T) -> DP(zcut=2.0)
  Path B: PARSI -> FILTER95(zcut1=1.0) -> pipe96(top 2) -> DALICON(lfitz=T) -> DP(zcut=2.0)

Both paths run in forward (query=cd1) and reverse (query=cd2) directions.
"""

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from wolf.wolf import wolf_compare, setup_protein, load_protein, SpatialHashGrid, NEIBORCUTOFF
from dp.dp import run_dp, DCCPResult, _format_ranges
from dalicon.dalicon import run_dalicon
from parsi.parsi import dowork_parsi, init_parsi
from filter95.filter95 import run_filter95

PIPENBEST = 2  # max entries per (cd1,cd2) key from FILTER95->DALICON


# ---------------------------------------------------------------------------
# Format converters
# ---------------------------------------------------------------------------

def dccp_to_dalicon_lines(dccp_results):
    """Convert DP DCCPResult list to DALICON input format string.

    Replaces dccp2dalicon.pl. Groups results by cd1 with END delimiters.
    """
    lines = []
    old_cd1 = None

    for r in dccp_results:
        lines.append('')  # blank line before each entry
        if r.cd1 != old_cd1:
            old_cd1 = r.cd1
            lines.append('END')
            lines.append(r.cd1)
        lines.append(f'{r.cd2}*')
        lines.append(str(r.nblock))
        # cd1 ranges in Fortran 8(i4,2x,i4) format
        lines.extend(_format_ranges(r.l1, r.r1, r.nblock))
        # cd2 ranges
        lines.extend(_format_ranges(r.l2, r.r2, r.nblock))

    lines.append('')
    lines.append('END')
    return '\n'.join(lines) + '\n'


def filter95_to_dalicon_lines(filter95_lines, zcut=1.0, nbest=PIPENBEST):
    """Convert FILTER95 output lines to DALICON input format string.

    Replaces: sort -nr | uniq | pipe96-free.pl zcut nbest
    """
    # Parse lines
    entries = []
    seen_lines = set()
    for line in filter95_lines:
        stripped = line.strip()
        if not stripped:
            continue
        # Dedup (uniq on sorted input = remove exact duplicate lines)
        if stripped in seen_lines:
            continue
        seen_lines.add(stripped)

        # Remove -99 placeholders (matching Perl s/\-99//g)
        cleaned = stripped.replace('-99', '')
        parts = cleaned.split()
        if len(parts) < 5:
            continue
        try:
            zscore = float(parts[0])
        except ValueError:
            continue

        key = parts[1]  # cd1cd2 (10 chars)
        # parts[2]=idom, parts[3]=score
        nseg = int(parts[4])
        range_vals = [int(x) for x in parts[5:5 + 4 * nseg]]
        entries.append((zscore, key, nseg, range_vals))

    # Sort by z-score descending (matching sort -nr)
    entries.sort(key=lambda x: -x[0])

    # Group by key, keep top nbest per key
    key_counts = {}
    output_lines = []
    old_cd1 = '?????'

    for zscore, key, nseg, range_vals in entries:
        key_counts[key] = key_counts.get(key, 0) + 1
        if key_counts[key] > nbest:
            continue
        if zscore < zcut:
            continue

        cd1 = key[:5].strip()
        cd2 = key[5:].strip()

        if cd1 != old_cd1:
            if old_cd1 != '?????':
                output_lines.append('END')
            output_lines.append(cd1)
            old_cd1 = cd1

        n = nseg
        output_lines.append(f'{cd2}*')
        output_lines.append(str(n))
        output_lines.append(' '.join(str(v) for v in range_vals[:2 * n]))
        output_lines.append(' '.join(str(v) for v in range_vals[2 * n:4 * n]))

    output_lines.append('END')
    output_lines.append('END')
    return '\n'.join(output_lines) + '\n'


# ---------------------------------------------------------------------------
# Pipeline paths
# ---------------------------------------------------------------------------

def run_wolf_path(query, targets, dat_dir):
    """Wolf path: WOLF -> DP(zcut=2.0) -> DALICON(lfitz=T) -> DP(zcut=2.0).

    Args:
        query: structure code to build grid for
        targets: list of structure codes to compare against
        dat_dir: path to .dat files

    Returns:
        list of DCCPResult from the final DP
    """
    dat_dir = Path(dat_dir)

    # Step 1: WOLF
    cd1_data = setup_protein(dat_dir / f'{query}.dat')
    prot1 = cd1_data[0]

    wolf_results = []
    if prot1.nseg > 2:
        grid = SpatialHashGrid()
        load_protein(grid, prot1.nseg, cd1_data[1], cd1_data[2],
                     cd1_data[3], NEIBORCUTOFF)
        for target in targets:
            result = wolf_compare(query, target, str(dat_dir),
                                  grid=grid, cd1_data=cd1_data)
            if result is not None:
                wolf_results.append(result)

    if not wolf_results:
        return []

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Step 2: DP on WOLF output
        wolfitz_path = tmpdir / 'wolfitz.txt'
        with open(wolfitz_path, 'w') as f:
            for r in wolf_results:
                f.write(r.to_wolfitz_line() + '\n')

        dp_results = run_dp(str(wolfitz_path), str(dat_dir), zcut=2.0)
        if not dp_results:
            return []

        # Step 3: DCCP -> DALICON input
        dalicon_text = dccp_to_dalicon_lines(dp_results)
        dalicon_path = tmpdir / 'dalicon_input.txt'
        with open(dalicon_path, 'w') as f:
            f.write(dalicon_text)

        # Step 4: DALICON (lfitz=True by default)
        dalicon_output = run_dalicon(str(dalicon_path), str(dat_dir))
        if not dalicon_output:
            return []

        # Step 5: Final DP
        wolfitz2_path = tmpdir / 'wolfitz2.txt'
        with open(wolfitz2_path, 'w') as f:
            for line in dalicon_output:
                f.write(line + '\n')

        final_results = run_dp(str(wolfitz2_path), str(dat_dir), zcut=2.0)

    return final_results


def run_parsi_path(query, targets, dat_dir):
    """Parsi path: PARSI -> FILTER95 -> pipe96 -> DALICON -> DP.

    Args:
        query: structure code (search tree built for this)
        targets: list of structure codes to compare against
        dat_dir: path to .dat files

    Returns:
        list of DCCPResult from the final DP
    """
    dat_dir = Path(dat_dir)
    dat_dir_slash = str(dat_dir) + '/'  # filter95 needs trailing slash

    # Step 1: PARSI
    init_parsi()
    all_refine_lines = []
    cd1_cache = {}
    for target in targets:
        lines = dowork_parsi(query, target, str(dat_dir),
                             lfirstonly=True, cd1_cache=cd1_cache)
        all_refine_lines.extend(lines)

    if not all_refine_lines:
        return []

    # Step 2: FILTER95
    filter95_output = run_filter95(all_refine_lines, dat_dir_slash, zcut1=1.0)
    if not filter95_output:
        return []

    # Step 3: pipe96 -> DALICON input
    dalicon_text = filter95_to_dalicon_lines(filter95_output, zcut=1.0,
                                             nbest=PIPENBEST)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Step 4: DALICON
        dalicon_path = tmpdir / 'dalicon_input.txt'
        with open(dalicon_path, 'w') as f:
            f.write(dalicon_text)

        dalicon_output = run_dalicon(str(dalicon_path), str(dat_dir))
        if not dalicon_output:
            return []

        # Step 5: Final DP
        wolfitz_path = tmpdir / 'wolfitz.txt'
        with open(wolfitz_path, 'w') as f:
            for line in dalicon_output:
                f.write(line + '\n')

        final_results = run_dp(str(wolfitz_path), str(dat_dir), zcut=2.0)

    return final_results


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

def compare_pair(cd1, cd2, dat_dir):
    """Full comparison of two structures.

    Runs wolf and parsi paths in both forward (query=cd1) and reverse
    (query=cd2) directions, accumulates all DCCP results.

    Args:
        cd1: first structure code
        cd2: second structure code
        dat_dir: path to .dat files

    Returns:
        list of DCCPResult (self-comparisons filtered out)
    """
    all_results = []
    targets = [cd1, cd2]

    # Forward direction: query=cd1
    all_results.extend(run_wolf_path(cd1, targets, dat_dir))
    all_results.extend(run_parsi_path(cd1, targets, dat_dir))

    # Reverse direction: query=cd2
    all_results.extend(run_wolf_path(cd2, targets, dat_dir))
    all_results.extend(run_parsi_path(cd2, targets, dat_dir))

    # Filter out self-comparisons
    all_results = [r for r in all_results if r.cd1 != r.cd2]

    return all_results
