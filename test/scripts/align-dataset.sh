#!/usr/bin/env bash

set -euo pipefail

#############################################################
#  BWAmeth Alignment and Sorted BAM Generation Script
#
#  Usage:
#    ./align-dataset.sh <REF_FASTA> <FASTQ_DIR> <OUTPUT_DIR> [THREADS]
#
#  This script:
#   - Aligns all paired-end FASTQ files in <FASTQ_DIR> (including subfolders!) to <REF_FASTA>
#   - Handles both _R1/_R2 and _1/_2 naming schemes automatically
#   - Sorts the alignments (coordinate sort)
#   - Outputs sorted BAMs and their index files to <OUTPUT_DIR>
#
#  Arguments:
#    REF_FASTA:   Reference fasta (indexed with BWAmeth) (.fa)
#    FASTQ_DIR:   Directory containing FASTQs (can have subfolders)
#    OUTPUT_DIR:  Output directory for BAM and BAI files
#    THREADS:     (optional) Number of threads for BWA/samtools [default: 8]
#############################################################

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <REF_FASTA> <FASTQ_DIR> <OUTPUT_DIR> [THREADS]"
  echo "  Aligns and sorts paired-end FASTQ files with BWAmeth (supports subfolders and both _R1/_1 naming)."
  exit 1
fi

REF_FASTA=$(readlink -f "$1")
FASTQ_DIR=$(readlink -f "$2")
OUTPUT_DIR=$(readlink -f "$3")
THREADS="${4:-8}"

mkdir -p "$OUTPUT_DIR"

# Find all possible R1 fastqs in all subdirs
R1_FILES=$(find "$FASTQ_DIR" -type f \( -name '*_R1*.fastq.gz' -o -name '*_1*.fastq.gz' \) | sort)
[[ -z "$R1_FILES" ]] && { echo "(!) No FASTQ R1 files found in $FASTQ_DIR"; exit 1; }

found_any=0
while IFS= read -r fq1; do
  # Determine pair naming scheme and get fq2
  if [[ "$fq1" =~ _R1 ]]; then
    fq2="${fq1/_R1/_R2}"
    sample=$(basename "$fq1" | sed 's/_R1.*\.fastq\.gz//')
  elif [[ "$fq1" =~ _1 ]]; then
    fq2="${fq1/_1/_2}"
    sample=$(basename "$fq1" | sed 's/_1.*\.fastq\.gz//')
  else
    echo "(!) Unrecognized FASTQ pair format: $fq1"
    continue
  fi

  if [[ ! -f "$fq2" ]]; then
    echo "(!) Missing pair for $sample ($fq2) — skipping"
    continue
  fi

  found_any=1
  echo "[*] Aligning and sorting $sample with BWAmeth..."
  bam_out="${OUTPUT_DIR}/${sample}.sorted.bam"

  bwameth.py --reference "$REF_FASTA" --threads "$THREADS" "$fq1" "$fq2" | \
    samtools sort -@ "$THREADS" -o "$bam_out" -

  samtools index "$bam_out"
  echo "[+] Done: $bam_out"
done <<< "$R1_FILES"

if [[ $found_any -eq 0 ]]; then
  echo "(!) No valid FASTQ pairs processed!"
  exit 1
fi

echo "========================================"
echo "All FASTQ pairs aligned, sorted, and indexed!"
