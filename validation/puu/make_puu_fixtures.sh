#!/bin/bash
# Fortran reference outputs for the puu port: for each structure run DaliLite's own
# dsspcmbi -> puu -> puutos exactly as mpidali.pm prepare_brkfile does, and keep
#   <id>/<id>.dssp   DSSP input to puu (fed to the Rust port, isolating it from DSSP)
#   <id>/units.puu   puu output tree (what the Rust port must reproduce)
#   <id>/<id>.dat    puutos output
# Usage: make_puu_fixtures.sh <pdb[.gz]> <4-char code> <chain> [...triplets]
set -e
DALI=/home/rschaeff/src/Dali_v5/DaliLite.v5/bin
OUT=$(cd "$(dirname "$0")" && pwd)
while [ $# -ge 3 ]; do
  pdb=$1; code=$2; chain=$3; shift 3
  id=$code$chain
  w=$(mktemp -d /tmp/puufx.XXXX)          # short path: DaliLite has an 80-char limit
  case $pdb in *.gz) gzip -cd "$pdb" > $w/x.pdb ;; *) cp "$pdb" $w/x.pdb ;; esac
  ( cd $w
    $DALI/dsspcmbi x.pdb $code.dssp > /dev/null 2>&1
    printf '#lhetatm\nTRUE\n#printdepth\n-99\n#lsplit\nTRUE\n' > puu.default
    $DALI/puu $id x.pdb $code.dssp > puu.log 2>&1
    $DALI/puutos units.puu $code.dssp ./ > puutos.log 2>&1 )
  mkdir -p $OUT/$id
  cp $w/$code.dssp $w/units.puu $OUT/$id/
  cp $w/$id.dat $OUT/$id/ 2>/dev/null || echo "no .dat for $id"
  rm -rf $w
done
