#!/usr/bin/env bash
set -euo pipefail

#############################################################################
# Region-targeted FASTQ files
#
# USAGE
#   bash create-test-fastq-files.sh <BAM_DIR> <TARGETS.bed> <OUT_DIR>
#
#############################################################################

[[ $# -lt 3 ]] && {
  echo "Usage: $0 <BAM_DIR> <TARGETS.bed> <SAMPLESHEET.csv> <OUT_DIR>" >&2
  exit 1
}

BAM_DIR=$(readlink -f "$1")
TARGETS=$(readlink -f "$2")
SAMPLES=$(readlink -f "$3")
OUT_DIR=$(readlink -f "$4")

THREADS="${5:-8}"
PAD="${6:-100}"

mkdir -p "$OUT_DIR"

shopt -s nullglob

echo "sample_id,read1,read2" > $SAMPLES

for BAM in "$BAM_DIR"/*.bam; do
    SAMPLE=$(basename "$BAM" | cut -d. -f1)
    FILTERED_BAM="${OUT_DIR}/${SAMPLE}.namesort.bam"

    echo -e "\n========= $SAMPLE ========="

    ## select reads overlapped bed file regions and then sort the result
    samtools view -b -L ${TARGETS} ${BAM} | \
        samtools sort -n -@ $THREADS -o ${FILTERED_BAM} -

    L1R1=${OUT_DIR}/${SAMPLE}_L001_1.fastq.gz
    L1R2=${OUT_DIR}/${SAMPLE}_L001_2.fastq.gz
    L2R1=${OUT_DIR}/${SAMPLE}_L002_1.fastq.gz
    L2R2=${OUT_DIR}/${SAMPLE}_L002_2.fastq.gz
    
    ## extract reads to generate paired-end fastq files
    samtools fastq -@ $THREADS -1 $L1R1 -2 $L1R2 \
        -0 /dev/null -s /dev/null ${FILTERED_BAM}

    cp $L1R1 $L2R1
    cp $L1R2 $L2R2

    L1R1=data/raw/`basename $L1R1`
    L1R2=data/raw/`basename $L1R2`
    L2R1=data/raw/`basename $L2R1`
    L2R2=data/raw/`basename $L2R2`
    
    echo "${SAMPLE},${L1R1},${L1R2}" >> $SAMPLES
    echo "${SAMPLE},${L2R1},${L2R2}" >> $SAMPLES
    
    rm ${FILTERED_BAM}
done

shopt -u nullglob

echo -e "Fastq files generated in $OUT_DIR"
