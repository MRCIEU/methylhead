#!/bin/bash

set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <files.csv> <sites.csv>"
  exit 1
fi

FILES="$1"
REGIONS="$2"

if [[ ! -f "$FILES" ]]; then
  echo "Error: File '$FILES' not found."
  exit 1
fi

MIN_COV=10
MIN_SAMPLES=$(($(wc -l < $FILES) / 2))

SAMPLE_FILES=$(mktemp)
trap 'rm -f "$SAMPLE_FILES"' EXIT

echo "collect-entropy-regions.sh---------------------------"
echo "$FILES $MIN_COV $MIN_SAMPLES $SAMPLE_FILES $REGIONS"

cut -d , -f 2 $FILES | tail -n +2 > $SAMPLE_FILES
csvtk concat -X $SAMPLE_FILES \
    | csvtk filter -f "depth>=$MIN_COV" \
    | csvtk cut -f 1,2,3 \
    | csvtk freq -f 1,2,3 \
    | csvtk filter -f "frequency>=$MIN_SAMPLES" \
    | (head -n 1 && tail -n +2 | sort -t, -k1,1 -k2,2n) \
    > $REGIONS
