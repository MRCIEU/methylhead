#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Region-targeted BAM ➜ FASTQ (Apptainer only)
#
# USAGE
#   bash create-test-fastq-files.sh <BAM_DIR> <TARGETS.bed> <OUT_DIR>
#                                [THREADS=8] [PAD=100]
#
###############################################################################

[[ $# -lt 4 ]] && {
  echo "Usage: $0 <BAM_DIR> <TARGETS.bed> <OUT_DIR> [THREADS] [PAD]" >&2
  exit 1
}

BAM_DIR=$(readlink -f "$1")
ORIG_BED=$(readlink -f "$2")
OUT_DIR=$(readlink -f "$3")

mkdir -p "$OUT_DIR"

THREADS="${5:-8}"
PAD="${6:-100}"

##################################################
## container config
##################################################
SIF_NAME="methylhead-pipeline_wgbs_analysis.sif"
SIF_ORAS="oras://docker.io/onuroztornaci/methylhead-pipeline:wgbs_analysis"

rm -f "$SIF_NAME"
apptainer pull "$SIF_NAME" "$SIF_ORAS"

###############################################################################
# 
###############################################################################
shopt -s nullglob

for BAM in "$BAM_DIR"/*.bam; do
    SAMPLE=$(basename "$BAM" | cut -d. -f1)
    FILTERED_BAM="${OUT_DIR}/${SAMPLE}.namesort.bam"

    echo -e "\n========= $SAMPLE ========="

    apptainer exec \
	--bind ${ORIG_BED}:${ORIG_BED} \
	--bind ${BAM_DIR}:${BAM_DIR} \
	--bind ${OUT_DIR}:${OUT_DIR} \
	"$SIF_NAME" \
	bash -euo pipefail -c "
	samtools view -b -L ${ORIG_BED} ${BAM_DIR} | \
        samtools sort -n -@ $THREADS -o ${FILTERED_BAM} -
        "

    apptainer exec \
	--bind ${OUT_DIR}:${OUT_DIR} \
	"$SIF_NAME" \
	bash -euo pipefile -c "
        samtools fastq -@ $THREADS \
        -1 ${OUT_DIR}/${SAMPLE}_1.fastq.gz \
        -2 ${OUT_DIR}/${SAMPLE}_2.fastq.gz \
        -0 /dev/null -s /dev/null ${FILTERED_BAM}
        "
    
    rm ${FILTERED_BAM}
done

shopt -u nullglob

do_cleanup() {
  echo -e "\n[*] Final cleanup…"
  rm -f "$SIF_NAME"
}

trap do_cleanup EXIT
echo -e "Fastq files generated in $OUT_DIR"
