#!/bin/bash
# Build and run the OPAL 1994 data/MC validation on grendel01.
#   ./run_grendel01.sh [output-dir] [extra CompareDataMC options...]
# e.g. ./run_grendel01.sh output --maxEvents 100000
set -euo pipefail

DIR=$(cd "$(dirname "$0")" && pwd)
MERGED=/raid5/data/yjlee/OPAL/converted/1994/merged
OUT=$DIR/output
# First argument is the output directory unless it is an option.
if [ $# -gt 0 ] && [ "${1#--}" = "$1" ]; then
  OUT=$1
  shift
fi

mkdir -p "$OUT"
make -C "$DIR"
"$DIR/CompareDataMC" \
  --data "$MERGED/OPAL_1994_data.root" \
  --mc "$MERGED/OPAL_1994_mc_jt74mh.root" \
  --output "$OUT/OPAL_1994_DataMC_validation.root" \
  --plotDir "$OUT" \
  "$@"
