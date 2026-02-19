#!/bin/bash
#SBATCH --job-name=dali_accuracy
#SBATCH --output=/home/rschaeff/dev/dali_cl/validation/batch01/accuracy_%j.log
#SBATCH --error=/home/rschaeff/dev/dali_cl/validation/batch01/accuracy_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=All

# Rust DALI accuracy test: batch_01 (999 proteins)
# Stage 1: Copy ECOD70 PDBs to /scratch (local I/O)
# Stage 2: Pre-generate .dat files from local PDBs
# Stage 3: Run accuracy comparison with cached .dat files

source ~/.bashrc
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libgomp.so.1

MANIFEST=/home/rschaeff/dev/dali_cl/validation/batch01/manifest.json
SCRIPT=/home/rschaeff/dev/dali_cl/test_batch_comparison.py
NFS_ECOD70=/home/rschaeff/data/dpam_reference/ecod_data/ECOD70
LOCAL_ECOD70=/scratch/dali_ecod70_$$
DAT_DIR=/scratch/dali_dat_cache_$$

echo "Starting accuracy test: $(date)"
echo "Host: $(hostname)"
echo "Python: $(which python3)"
echo "Manifest: $MANIFEST"
echo ""

# Stage 1: Copy ECOD70 PDBs to local scratch
echo "Stage 1: Copying ECOD70 to $LOCAL_ECOD70..."
t0=$(date +%s)
mkdir -p "$LOCAL_ECOD70"
rsync -a "$NFS_ECOD70/" "$LOCAL_ECOD70/"
t1=$(date +%s)
echo "  Copied in $((t1-t0))s ($(ls $LOCAL_ECOD70/*.pdb | wc -l) PDB files)"
echo ""

# Stage 2: Pre-generate .dat files
echo "Stage 2: Pre-generating .dat files to $DAT_DIR..."
mkdir -p "$DAT_DIR"
python3 -u -c "
import json, sys, time
sys.path.insert(0, '/home/rschaeff/dev/dpam_c2')
import dali as dali_rust
from pathlib import Path

manifest = json.load(open('$MANIFEST'))
ecod70 = Path('$LOCAL_ECOD70')
dat_dir = Path('$DAT_DIR')

# Collect unique templates
templates = set()
for p in manifest:
    with open(p['hits4dali']) as f:
        for line in f:
            t = line.strip()
            if t:
                templates.add(t)
print(f'Unique templates: {len(templates)}')

t0 = time.time()
n_ok, n_fail = 0, 0
for i, code in enumerate(sorted(templates)):
    dat_path = dat_dir / f'{code}.dat'
    if dat_path.exists():
        n_ok += 1
        continue
    pdb_path = ecod70 / f'{code}.pdb'
    if not pdb_path.exists():
        n_fail += 1
        continue
    try:
        protein = dali_rust.import_pdb(str(pdb_path), '', code)
        dali_rust.write_dat(protein, str(dat_path))
        n_ok += 1
    except Exception as e:
        n_fail += 1
    if (i+1) % 5000 == 0:
        elapsed = time.time() - t0
        print(f'  {i+1}/{len(templates)} ({n_ok} ok, {n_fail} fail) {elapsed:.0f}s')

elapsed = time.time() - t0
print(f'Done: {n_ok} ok, {n_fail} fail in {elapsed:.0f}s')
"
echo ""

# Stage 3: Run accuracy test
echo "Stage 3: Running accuracy comparison..."
echo ""
python3 -u "$SCRIPT" "$MANIFEST" --dat-dir "$DAT_DIR"

echo ""
echo "Completed: $(date)"

# Cleanup local scratch
rm -rf "$LOCAL_ECOD70" "$DAT_DIR"
