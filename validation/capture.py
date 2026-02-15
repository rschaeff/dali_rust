#!/usr/bin/env python3
"""
Capture ground truth from DaliLite.v5 at each pipeline stage.

Runs serialcompare for each method independently, saving all inputs
and outputs in a structured fixtures directory. These fixtures serve
as the validation baseline for any reimplementation.

Usage:
    python capture.py [--dali-root /path/to/DaliLite.v5] [--dat-dir /path/to/DAT]
"""

import argparse
import json
import os
import shutil
import subprocess
import tempfile
import time
from itertools import product
from pathlib import Path


def get_structures(dat_dir):
    """List available structure IDs from .dat files."""
    return sorted([
        p.stem for p in Path(dat_dir).glob('*.dat')
    ])


def run_stage(cmd, stdin_file=None, work_dir=None, timeout=300):
    """Run a shell command, return stdout, stderr, elapsed time."""
    t0 = time.time()
    stdin_data = None
    if stdin_file and os.path.exists(stdin_file):
        with open(stdin_file) as f:
            stdin_data = f.read()

    result = subprocess.run(
        cmd, shell=True, cwd=work_dir,
        input=stdin_data, capture_output=True, text=True,
        timeout=timeout,
    )
    elapsed = time.time() - t0
    return result.stdout, result.stderr, result.returncode, elapsed


def collect_fort_files(work_dir):
    """Collect and concatenate fort.1[0-9][0-9] files, return content."""
    import glob
    pattern = os.path.join(work_dir, 'fort.1[0-9][0-9]')
    files = sorted(glob.glob(pattern))
    content = ''
    for f in files:
        with open(f) as fh:
            content += fh.read()
    return content, files


def cleanup_fort_files(work_dir):
    """Remove fort.* files from work directory."""
    import glob
    for f in glob.glob(os.path.join(work_dir, 'fort.*')):
        os.remove(f)


def capture_wolf(compare_exe, dat_dir, structures, fixtures_dir, work_dir):
    """Capture WOLF stage: SSE-based screening for all pairs."""
    stage_dir = Path(fixtures_dir) / 'wolf'
    stage_dir.mkdir(exist_ok=True)

    # Write list files
    list1_path = os.path.join(work_dir, 'list1')
    list2_path = os.path.join(work_dir, 'list2')
    with open(list1_path, 'w') as f:
        for s in structures:
            f.write(s + '\n')
    with open(list2_path, 'w') as f:
        for s in structures:
            f.write(s + '\n')

    cleanup_fort_files(work_dir)

    # Run WOLF
    cmd = f'{compare_exe} {dat_dir} {dat_dir} WOLF'
    stdout, stderr, rc, elapsed = run_stage(cmd, work_dir=work_dir)

    # Collect output
    output, fort_files = collect_fort_files(work_dir)

    # Save
    with open(stage_dir / 'wolf_output.txt', 'w') as f:
        f.write(output)
    with open(stage_dir / 'list1.txt', 'w') as f:
        for s in structures:
            f.write(s + '\n')
    with open(stage_dir / 'list2.txt', 'w') as f:
        for s in structures:
            f.write(s + '\n')
    with open(stage_dir / 'metadata.json', 'w') as f:
        json.dump({
            'command': cmd,
            'return_code': rc,
            'elapsed_seconds': elapsed,
            'n_structures': len(structures),
            'n_fort_files': len(fort_files),
            'structures': structures,
        }, f, indent=2)

    print(f'  WOLF: {len(output.splitlines())} output lines, {elapsed:.1f}s')
    cleanup_fort_files(work_dir)
    return stage_dir / 'wolf_output.txt'


def capture_dp(compare_exe, dat_dir, wolf_output, fixtures_dir, work_dir, label='wolf'):
    """Capture DP stage: Z-score computation from WOLFITZ input."""
    stage_dir = Path(fixtures_dir) / 'dp' / label
    stage_dir.mkdir(parents=True, exist_ok=True)

    cleanup_fort_files(work_dir)

    # Copy wolf_output to work_dir for stdin
    input_path = os.path.join(work_dir, 'dp_stdin')
    shutil.copy(wolf_output, input_path)

    cmd = f'{compare_exe} {dat_dir} {dat_dir} DP < {input_path}'
    stdout, stderr, rc, elapsed = run_stage(cmd, work_dir=work_dir)

    output, fort_files = collect_fort_files(work_dir)

    with open(stage_dir / 'dp_input.txt', 'w') as f:
        with open(wolf_output) as wf:
            f.write(wf.read())
    with open(stage_dir / 'dp_output.txt', 'w') as f:
        f.write(output)
    with open(stage_dir / 'metadata.json', 'w') as f:
        json.dump({
            'command': cmd,
            'return_code': rc,
            'elapsed_seconds': elapsed,
            'n_fort_files': len(fort_files),
        }, f, indent=2)

    print(f'  DP ({label}): {len(output.splitlines())} output lines, {elapsed:.1f}s')
    cleanup_fort_files(work_dir)
    return stage_dir / 'dp_output.txt'


def capture_dalicon(compare_exe, dat_dir, dalicon_input, fixtures_dir, work_dir,
                    label='wolf'):
    """Capture DALICON stage: genetic algorithm alignment refinement."""
    stage_dir = Path(fixtures_dir) / 'dalicon' / label
    stage_dir.mkdir(parents=True, exist_ok=True)

    cleanup_fort_files(work_dir)

    input_path = os.path.join(work_dir, 'dalicon_stdin')
    shutil.copy(dalicon_input, input_path)

    cmd = f'{compare_exe} {dat_dir} {dat_dir} DALICON T < {input_path}'
    stdout, stderr, rc, elapsed = run_stage(cmd, work_dir=work_dir)

    output, fort_files = collect_fort_files(work_dir)

    with open(stage_dir / 'dalicon_input.txt', 'w') as f:
        with open(dalicon_input) as di:
            f.write(di.read())
    with open(stage_dir / 'dalicon_output.txt', 'w') as f:
        f.write(output)
    with open(stage_dir / 'metadata.json', 'w') as f:
        json.dump({
            'command': cmd,
            'return_code': rc,
            'elapsed_seconds': elapsed,
            'n_fort_files': len(fort_files),
        }, f, indent=2)

    print(f'  DALICON ({label}): {len(output.splitlines())} output lines, {elapsed:.1f}s')
    cleanup_fort_files(work_dir)
    return stage_dir / 'dalicon_output.txt'


def capture_parsi(compare_exe, dat_dir, structures, fixtures_dir, work_dir):
    """Capture PARSI stage: exhaustive branch-and-bound alignment."""
    stage_dir = Path(fixtures_dir) / 'parsi'
    stage_dir.mkdir(exist_ok=True)

    # PARSI reads list1/list2 from work_dir
    list1_path = os.path.join(work_dir, 'list1')
    list2_path = os.path.join(work_dir, 'list2')
    with open(list1_path, 'w') as f:
        for s in structures:
            f.write(s + '\n')
    with open(list2_path, 'w') as f:
        for s in structures:
            f.write(s + '\n')

    cleanup_fort_files(work_dir)

    # PARSI also accepts pairs on stdin (PARSI12 mode), but the standard
    # mode reads list1/list2. Let's capture both modes.

    # Standard PARSI (all-vs-all from list files)
    cmd = f'{compare_exe} {dat_dir} {dat_dir} PARSI'
    stdout, stderr, rc, elapsed = run_stage(cmd, work_dir=work_dir, timeout=600)

    output, fort_files = collect_fort_files(work_dir)

    with open(stage_dir / 'parsi_output.txt', 'w') as f:
        f.write(output)
    with open(stage_dir / 'list1.txt', 'w') as f:
        for s in structures:
            f.write(s + '\n')
    with open(stage_dir / 'list2.txt', 'w') as f:
        for s in structures:
            f.write(s + '\n')
    with open(stage_dir / 'metadata.json', 'w') as f:
        json.dump({
            'command': cmd,
            'return_code': rc,
            'elapsed_seconds': elapsed,
            'n_structures': len(structures),
            'n_fort_files': len(fort_files),
        }, f, indent=2)

    print(f'  PARSI: {len(output.splitlines())} output lines, {elapsed:.1f}s')
    cleanup_fort_files(work_dir)
    return stage_dir / 'parsi_output.txt'


def capture_filter95(compare_exe, dat_dir, parsi_output, fixtures_dir, work_dir):
    """Capture FILTER95 stage: filter and refine PARSI output."""
    stage_dir = Path(fixtures_dir) / 'filter95'
    stage_dir.mkdir(exist_ok=True)

    cleanup_fort_files(work_dir)

    input_path = os.path.join(work_dir, 'filter_stdin')
    shutil.copy(parsi_output, input_path)

    cmd = f'{compare_exe} {dat_dir} {dat_dir} FILTER95 < {input_path}'
    stdout, stderr, rc, elapsed = run_stage(cmd, work_dir=work_dir)

    output, fort_files = collect_fort_files(work_dir)

    with open(stage_dir / 'filter95_input.txt', 'w') as f:
        with open(parsi_output) as pi:
            f.write(pi.read())
    with open(stage_dir / 'filter95_output.txt', 'w') as f:
        f.write(output)
    with open(stage_dir / 'stderr.txt', 'w') as f:
        f.write(stderr)
    with open(stage_dir / 'metadata.json', 'w') as f:
        json.dump({
            'command': cmd,
            'return_code': rc,
            'elapsed_seconds': elapsed,
            'n_fort_files': len(fort_files),
        }, f, indent=2)

    print(f'  FILTER95: {len(output.splitlines())} output lines, {elapsed:.1f}s')
    cleanup_fort_files(work_dir)
    return stage_dir / 'filter95_output.txt'


def capture_e2e_pairwise(dali_root, dat_dir, cd1, cd2, fixtures_dir, work_dir):
    """Capture end-to-end pairwise comparison result."""
    stage_dir = Path(fixtures_dir) / 'e2e' / f'{cd1}_vs_{cd2}'
    stage_dir.mkdir(parents=True, exist_ok=True)

    pairwise_dir = os.path.join(work_dir, f'e2e_{cd1}_{cd2}')
    os.makedirs(pairwise_dir, exist_ok=True)

    dali_pl = os.path.join(dali_root, 'bin', 'dali.pl')
    cmd = (f'perl {dali_pl} --cd1 {cd1} --cd2 {cd2}'
           f' --dat1 {dat_dir} --dat2 {dat_dir}'
           f' --outfmt "summary,alignments"')

    stdout, stderr, rc, elapsed = run_stage(cmd, work_dir=pairwise_dir, timeout=300)

    # Collect outputs
    for ext in ('dccp', 'txt', 'html'):
        outfile = os.path.join(pairwise_dir, f'{cd1}.{ext}')
        if os.path.exists(outfile):
            shutil.copy(outfile, stage_dir / f'{cd1}.{ext}')

    with open(stage_dir / 'stdout.txt', 'w') as f:
        f.write(stdout)
    with open(stage_dir / 'stderr.txt', 'w') as f:
        f.write(stderr)
    with open(stage_dir / 'metadata.json', 'w') as f:
        json.dump({
            'command': cmd,
            'return_code': rc,
            'elapsed_seconds': elapsed,
            'cd1': cd1,
            'cd2': cd2,
        }, f, indent=2)

    # Clean up lock file
    lock = os.path.join(pairwise_dir, 'dali.lock')
    if os.path.exists(lock):
        os.remove(lock)

    return stage_dir


def convert_dp_to_dalicon_input(dp_output, dalicon_input_path, dali_root):
    """Run dccp2dalicon.pl to convert DP DCCP output to DALICON input."""
    dccp2dalicon = os.path.join(dali_root, 'bin', 'dccp2dalicon.pl')
    cmd = f'perl {dccp2dalicon} < {dp_output} > {dalicon_input_path}'
    subprocess.run(cmd, shell=True, check=True)
    return dalicon_input_path


def convert_filter95_to_dalicon_input(filter95_output, dalicon_input_path, dali_root):
    """Run pipe96-free.pl to convert FILTER95 output to DALICON input."""
    pipe96 = os.path.join(dali_root, 'bin', 'pipe96-free.pl')
    cmd = (f'sort -nr {filter95_output} | uniq'
           f' | perl {pipe96} 1.0 2 > {dalicon_input_path}')
    subprocess.run(cmd, shell=True, check=True)
    return dalicon_input_path


def main():
    parser = argparse.ArgumentParser(description='Capture DaliLite ground truth')
    parser.add_argument('--dali-root',
                        default=os.path.expanduser('~/src/Dali_v5/DaliLite.v5'),
                        help='Path to DaliLite.v5 root')
    parser.add_argument('--dat-dir', default=None,
                        help='Path to DAT directory (default: dali-root/DAT)')
    parser.add_argument('--fixtures-dir',
                        default=os.path.join(os.path.dirname(__file__), 'fixtures'),
                        help='Output directory for fixtures')
    parser.add_argument('--structures', nargs='*', default=None,
                        help='Specific structures to test (default: all in DAT)')
    parser.add_argument('--stages', nargs='*',
                        default=['wolf', 'dp', 'parsi', 'filter95', 'dalicon', 'e2e'],
                        help='Stages to capture')
    parser.add_argument('--e2e-pairs', nargs='*', default=None,
                        help='Specific pairs for e2e (e.g., 2nrmA:101mA)')
    args = parser.parse_args()

    dali_root = args.dali_root
    dat_dir = args.dat_dir or os.path.join(dali_root, 'DAT')
    fixtures_dir = args.fixtures_dir
    compare_exe = os.path.join(dali_root, 'bin', 'serialcompare')

    # Resolve to absolute paths
    dat_dir = os.path.abspath(dat_dir) + '/'
    fixtures_dir = os.path.abspath(fixtures_dir)

    os.makedirs(fixtures_dir, exist_ok=True)

    # Get structure list
    structures = args.structures or get_structures(dat_dir)
    print(f'Structures: {structures}')
    print(f'DAT dir: {dat_dir}')
    print(f'Fixtures dir: {fixtures_dir}')
    print()

    # Copy .dat files to fixtures for reference
    struct_dir = Path(fixtures_dir) / 'structures'
    struct_dir.mkdir(exist_ok=True)
    for s in structures:
        src = os.path.join(dat_dir, f'{s}.dat')
        dst = struct_dir / f'{s}.dat'
        if not dst.exists():
            shutil.copy(src, dst)

    # Create a temp work directory
    work_dir = tempfile.mkdtemp(prefix='dali_capture_')
    print(f'Work directory: {work_dir}')
    print()

    # Copy list files to work_dir (WOLF/PARSI read from cwd)
    list1_path = os.path.join(work_dir, 'list1')
    list2_path = os.path.join(work_dir, 'list2')

    try:
        # Stage 1: WOLF
        if 'wolf' in args.stages:
            print('Stage 1: WOLF (SSE screening)')
            wolf_output = capture_wolf(
                compare_exe, dat_dir, structures, fixtures_dir, work_dir)
        else:
            wolf_output = Path(fixtures_dir) / 'wolf' / 'wolf_output.txt'

        # Stage 2: DP on WOLF output
        if 'dp' in args.stages and wolf_output.exists():
            print('Stage 2: DP (Z-scores from WOLF)')
            dp_output = capture_dp(
                compare_exe, dat_dir, wolf_output, fixtures_dir, work_dir,
                label='wolf')

            # Convert to DALICON input
            dalicon_input_wolf = os.path.join(work_dir, 'dalicon_input_wolf')
            convert_dp_to_dalicon_input(dp_output, dalicon_input_wolf, dali_root)
            dalicon_wolf_dir = Path(fixtures_dir) / 'dalicon' / 'wolf'
            dalicon_wolf_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy(dalicon_input_wolf,
                        dalicon_wolf_dir / 'dalicon_input_from_dp.txt')
        else:
            dp_output = None

        # Stage 3: DALICON on WOLF->DP output
        if 'dalicon' in args.stages and dp_output and os.path.exists(
                os.path.join(work_dir, 'dalicon_input_wolf')):
            print('Stage 3: DALICON (alignment refinement from WOLF)')
            dalicon_input_wolf = os.path.join(work_dir, 'dalicon_input_wolf')
            dalicon_output = capture_dalicon(
                compare_exe, dat_dir, dalicon_input_wolf, fixtures_dir, work_dir,
                label='wolf')

            # DP on DALICON output
            print('Stage 3b: DP (Z-scores from DALICON)')
            dp_dalicon_output = capture_dp(
                compare_exe, dat_dir, dalicon_output, fixtures_dir, work_dir,
                label='dalicon')

        # Stage 4: PARSI
        if 'parsi' in args.stages:
            print('Stage 4: PARSI (exhaustive alignment)')
            parsi_output = capture_parsi(
                compare_exe, dat_dir, structures, fixtures_dir, work_dir)
        else:
            parsi_output = Path(fixtures_dir) / 'parsi' / 'parsi_output.txt'

        # Stage 5: FILTER95
        if 'filter95' in args.stages and parsi_output.exists():
            print('Stage 5: FILTER95 (filter and refine)')
            filter95_output = capture_filter95(
                compare_exe, dat_dir, parsi_output, fixtures_dir, work_dir)

            # Convert to DALICON input
            dalicon_input_parsi = os.path.join(work_dir, 'dalicon_input_parsi')
            convert_filter95_to_dalicon_input(
                filter95_output, dalicon_input_parsi, dali_root)

            parsi_dalicon_dir = Path(fixtures_dir) / 'dalicon' / 'parsi'
            parsi_dalicon_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy(dalicon_input_parsi,
                        parsi_dalicon_dir / 'dalicon_input_from_parsi.txt')

            # DALICON on PARSI path
            if 'dalicon' in args.stages:
                print('Stage 5b: DALICON (alignment refinement from PARSI)')
                capture_dalicon(
                    compare_exe, dat_dir, dalicon_input_parsi, fixtures_dir,
                    work_dir, label='parsi')

        # Stage 6: End-to-end pairwise comparisons
        if 'e2e' in args.stages:
            print('Stage 6: End-to-end pairwise comparisons')
            if args.e2e_pairs:
                pairs = [p.split(':') for p in args.e2e_pairs]
            else:
                # Default: representative pairs covering different scenarios
                pairs = []
                for i, s1 in enumerate(structures):
                    for s2 in structures[i:]:
                        pairs.append((s1, s2))

            for cd1, cd2 in pairs:
                print(f'  E2E: {cd1} vs {cd2}')
                try:
                    capture_e2e_pairwise(
                        dali_root, dat_dir, cd1, cd2, fixtures_dir, work_dir)
                except Exception as e:
                    print(f'    FAILED: {e}')

    finally:
        print(f'\nWork directory preserved at: {work_dir}')
        print(f'Fixtures saved to: {fixtures_dir}')

    # Write summary
    summary = {
        'dali_root': dali_root,
        'dat_dir': dat_dir,
        'structures': structures,
        'stages_captured': args.stages,
        'compare_exe': compare_exe,
    }
    with open(os.path.join(fixtures_dir, 'capture_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)

    print('\nCapture complete.')


if __name__ == '__main__':
    main()
