#!/usr/bin/env bash
set -euo pipefail

##################################################
# Select a subset of regions from lung cancer panel 
# as targets for this small dataset
#
# USAGE
# bash select-test-panel.sh <BAM_DIR> <PANEL_CSV> \
#       <CELL_BED> <OUT_CSV> <OUT_BED>
##################################################

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 <BAM_DIR> <PANEL_CSV> <CELL_BED> <OUT_CSV> <OUT_BED>"
  exit 1
fi

BAM_DIR=$(readlink -f "$1")
PANEL_CSV=$(readlink -f "$2")
CELL_BED=$(readlink -f "$3")
OUT_CSV=$(readlink -f "$4")
OUT_BED=$(readlink -f "$5")

##################################################
## container config
##################################################
SIF_NAME="methylhead-pipeline_wgbs_analysis.sif"
SIF_ORAS="oras://docker.io/onuroztornaci/methylhead-pipeline:wgbs_analysis"

rm -f "$SIF_NAME"
apptainer pull "$SIF_NAME" "$SIF_ORAS"

##################################################
## count reads in bam files overlapping bed file regions
##################################################
apptainer exec \
    --bind "$CELL_BED":"$CELL_BED" \
    --bind "$BAM_DIR":"$BAM_DIR" \
    "$SIF_NAME" \
    samtools bedcov $CELL_BED $BAM_DIR/*.bam > coverage.bed

##################################################
## create panel
##################################################
Rscript scripts/select-test-panel.r $PANEL_CSV coverage.bed $OUT_CSV $OUT_BED

echo "Panel files generated: $OUT_CSV and $OUT_BED"

##################################################
## cleanup
##################################################
do_cleanup() {
  echo -e "\n[*] Final cleanup…"
  rm -f coverage.bed
  rm -f "$SIF_NAME" 
}

trap do_cleanup EXIT

