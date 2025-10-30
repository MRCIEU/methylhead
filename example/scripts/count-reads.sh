#!/usr/bin/env bash
set -euo pipefail

##################################################
# Calculate numbers of reads in bed file regions
#
# USAGE
# bash scripts/count-reads.sh <BAM_DIR> <IN_BED> <COUNTS_BED> 
##################################################


if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <BAM_DIR> <IN_BED> <COUNTS_BED>"
  exit 1
fi

BAM_DIR=$(readlink -f "$1")
IN_BED=$(readlink -f "$2")
COUNTS_BED=$(readlink -f "$3")

## count reads in bam files overlapping bed file regions
samtools bedcov $IN_BED $BAM_DIR/*.bam > $COUNTS_BED

