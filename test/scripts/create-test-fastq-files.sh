#!/usr/bin/env bash
set -euo pipefail

#############################################################################
# Region-targeted BAM ➜ FASTQ
#
# USAGE
#   bash create-test-fastq-files.sh <BAM_DIR> <TARGETS.bed> <OUT_DIR>
#                                [THREADS=8] [PAD=100]
#
#############################################################################

[[ $# -lt 3 ]] && {
  echo "Usage: $0 <BAM_DIR> <TARGETS.bed> <OUT_DIR> [THREADS] [PAD]" >&2
  exit 1
}

BAM_DIR=$(readlink -f "$1")
ORIG_BED=$(readlink -f "$2")
OUT_DIR=$(readlink -f "$3")

THREADS="${5:-8}"
PAD="${6:-100}"

mkdir -p "$OUT_DIR"

shopt -s nullglob

for BAM in "$BAM_DIR"/*.bam; do
    SAMPLE=$(basename "$BAM" | cut -d. -f1)
    FILTERED_BAM="${OUT_DIR}/${SAMPLE}.namesort.bam"

    echo -e "\n========= $SAMPLE ========="

    ## select reads overlapped bed file regions and then sort the result
    samtools view -b -L ${ORIG_BED} ${BAM} | \
        samtools sort -n -@ $THREADS -o ${FILTERED_BAM} -

    ## extract reads to generate paired-end fastq files
    samtools fastq -@ $THREADS \
        -1 ${OUT_DIR}/${SAMPLE}_1.fastq.gz \
        -2 ${OUT_DIR}/${SAMPLE}_2.fastq.gz \
        -0 /dev/null -s /dev/null ${FILTERED_BAM}
    
    rm ${FILTERED_BAM}
done

shopt -u nullglob

echo -e "Fastq files generated in $OUT_DIR"
