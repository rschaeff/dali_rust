#!/bin/bash
# Capture DALICON ground truth for the expanded corpus
# Pipeline: WOLF → DP → DALICON → DP
#
# Uses the original DaliLite Fortran binaries to generate reference output

set -e

DALI_BIN="/home/rschaeff/src/Dali_v5/DaliLite.v5/bin"
DAT_DIR="/home/rschaeff/dev/dali_cl/validation/corpus_expansion/dat"
WORK_DIR="/home/rschaeff/dev/dali_cl/validation/corpus_expansion/work"
OUT_DIR="/home/rschaeff/dev/dali_cl/validation/corpus_expansion/ground_truth"

mkdir -p "$OUT_DIR"/{wolf,dp_wolf,dalicon_wolf,dp_dalicon}
cd "$WORK_DIR"

# List of structures (chain codes matching .dat files)
STRUCTURES=(
    1a7sA 1aiwA 1a25A 1a12A 1f3uA 1a04A 1bbhA 1b3qA 1a17A
    1miwA 1aopA 1a8lA 1a6qA 1a06A 1b3oB 1a0cA 1a4iA 1bcoA
)

echo "=== Structures: ${#STRUCTURES[@]} ==="
echo "=== Total pairs: $((${#STRUCTURES[@]} * ${#STRUCTURES[@]})) ==="

# Step 1: Generate list files for serialcompare
# list1.txt and list2.txt contain the structure codes
printf '%s\n' "${STRUCTURES[@]}" > "$OUT_DIR/list1.txt"
printf '%s\n' "${STRUCTURES[@]}" > "$OUT_DIR/list2.txt"

# Step 2: Run WOLF stage
echo ""
echo "=== Step 2: WOLF ==="
rm -f fort.101
"$DALI_BIN/serialcompare" "$DAT_DIR/" "$DAT_DIR/" WOLF T \
    < <(paste -d'\n' "$OUT_DIR/list1.txt" "$OUT_DIR/list2.txt" | head -n $((${#STRUCTURES[@]} * 2))) \
    > "$OUT_DIR/wolf/wolf_stderr.txt" 2>&1 || true

# Actually, serialcompare reads list1 then list2 from stdin
# Let me check the format...
echo "Checking serialcompare input format..."

