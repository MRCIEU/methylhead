#!/usr/bin/env bash
set -euo pipefail

################################################################################
# create-test-genome.sh   (container-aware)
#
# Usage:
#   bash create-test-genome.sh <original_fasta> <regions_bed> <genome_fasta>
#
#   <original_fasta>  : Original genome FASTA
#   <regions_csv>     : Target regions BED
#   <genome_fasta>    : Test genome FASTA
################################################################################

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <original_fasta> <regions_bed> <genome_fasta>"
  exit 1
fi

ORIGINAL_FASTA="$(readlink -f "$1")"
REGIONS_BED="$(readlink -f "$2")"
GENOME_FASTA="$(readlink -f "$3")"

GENOME_DIR=`dirname "$GENOME_FASTA"`
mkdir -p $GENOME_DIR

###############################################################################
# container config
###############################################################################

SIF_NAME="methylhead-pipeline_bedtools.sif"
SIF_ORAS="docker://quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0"

###############################################################################
# generate fasta file for all target regions
###############################################################################

apptainer exec \
    --bind "$ORIGINAL_FASTA":"$ORIGINAL_FASTA" \
    --bind "$REGIONS_BED":"$REGIONS_BED" \
    --bind "$GENOME_DIR":"$GENOME_DIR" \
    "$SIF_NAME" \
    bedtools getfasta -fi "$ORIGINAL_FASTA" -bed "$REGIONS_BED" -fo "${GENOME_FASTA}.tmp"

###############################################################################
# merge all segments from the same chromosome into a single FASTA entry
###############################################################################

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

##################################################
## cleanup
##################################################
do_cleanup() {
  echo -e "\n[*] Final cleanup…"
  rm -f "${GENOME_FASTA}.tmp"
  rm -f "${SIF_NAME}"
}

trap do_cleanup EXIT


