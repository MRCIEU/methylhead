#!/usr/bin/env bash
set -euo pipefail

################################################################################
# create-test-genome.sh   (container-aware)
#
# Usage:
#   bash create-test-genome.sh <original_fasta> <regions_bed> <genome_fasta>
#
#   <original_fasta>  : Original genome FASTA
#   <regions_bed>     : Target regions BED
#   <genome_fasta>    : Test genome FASTA
################################################################################

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <original_fasta> <regions_bed> <genome_fasta>"
  exit 1
fi

echo "here"

ORIGINAL_FASTA="$1" 
REGIONS_BED="$2"
GENOME_FASTA="$3"

GENOME_DIR=`dirname "$GENOME_FASTA"`

echo "$GENOME_DIR"

mkdir -p $GENOME_DIR

echo "$REGIONS_BED"

# generate fasta file for all target regions
bedtools getfasta -fi "$ORIGINAL_FASTA" \
    -bed "$REGIONS_BED" -fo "${GENOME_FASTA}.tmp"

echo "$GENOME_FASTA"

# merge all segments from the same chromosome into a single FASTA entry
awk '
/^>/ {
  split($0, arr, "[:]");
  chr = substr(arr[1], 2);
  if (!seq[chr]) order[++n] = chr;
  curr_chr = chr;
  next;
}
{ seq[curr_chr] = seq[curr_chr] $0 }
END {
  for (i = 1; i <= n; i++) {
    print ">" order[i];
    s = seq[order[i]];
    for (j = 1; j <= length(s); j += 60)
      print substr(s, j, 60);
  }
}
' "${GENOME_FASTA}.tmp" > "${GENOME_FASTA}"

echo "Fasta file generated: $GENOME_FASTA"

## cleanup
do_cleanup() {
  echo -e "\n[*] Final cleanup…"
  rm -f "${GENOME_FASTA}.tmp"
}

trap do_cleanup EXIT


