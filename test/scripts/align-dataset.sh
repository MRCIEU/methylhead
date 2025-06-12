#!/usr/bin/env bash
set -euo pipefail

#############################################################################
# Align paired-end reads
#
# USAGE
#   bash scripts/align-dataset.sh <REF_FASTA> <RAW_DIR> <OUT_DIR> [THREADS] 
#############################################################################

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <REF_FASTA> <RAW_FASTQ_DIR> <OUT_DIR> [THREADS]" >&2
  exit 1
fi

REF_FASTA=$(readlink -f "$1")
RAW_DIR=$(readlink -f "$2")
OUT_DIR=$(readlink -f "$3")
THREADS="${4:-8}"

TRIM_DIR="${OUT_DIR}/trimmed"
BAM_DIR="${OUT_DIR}/bam"
mkdir -p "$TRIM_DIR" "$BAM_DIR"


## BWAmeth index (one-off)
if ! ls "${REF_FASTA}".*bwameth.c2t* &>/dev/null; then
    echo "[*] Building BWAmeth index"
    bwameth.py index "$REF_FASTA"
else
    echo "[·] BWAmeth index already present."
fi

## FASTQ ➜ Trim ➜ BWAmeth ➜ BAM

## make a list of first read files (wide variety of possible name schemes)
mapfile -t RAW_R1_LIST < <(find "$RAW_DIR" -type f \( \
      -name '*_R1*.fastq.gz' -o -name '*_R1*.fq.gz' \
      -o -name '*_1*.fastq.gz'  -o -name '*_1*.fq.gz' \
    \) | sort)

## quit if there are no such files
[[ ${#RAW_R1_LIST[@]} -eq 0 ]] && { echo "(!) No FASTQs found"; exit 1; }

## for each first read file
for RAW_R1 in "${RAW_R1_LIST[@]}"; do
    ## identify the paired read file
    if [[ "$RAW_R1" =~ _R1 ]]; then
	RAW_R2="${RAW_R1/_R1/_R2}"; 
	SAMPLE=$(basename "$RAW_R1" | sed -E 's/_R1[^.]*\.f(ast)?q\.gz//')
    else
	RAW_R2="${RAW_R1/_1/_2}";   
	SAMPLE=$(basename "$RAW_R1" | sed -E 's/_1[^.]*\.f(ast)?q\.gz//')
    fi
    [[ -f "$RAW_R2" ]] || { 
	echo "(!) R2 missing for $SAMPLE, skipping."; continue; 
    }

    ## trim read pairs
    TRIMMED_R1="${TRIM_DIR}/${SAMPLE}_val_1.fq.gz"
    TRIMMED_R2="${TRIM_DIR}/${SAMPLE}_val_2.fq.gz"
    if [[ -s "$TRIMMED_R1" && -s "$TRIMMED_R2" ]]; then
	echo "[·] $SAMPLE already trimmed."
    else
	echo "[*] Trimming $SAMPLE"
	trim_galore --paired --gzip --cores "$THREADS" \
            --basename "$SAMPLE" --output_dir "$TRIM_DIR" \
	    "$RAW_R1" "$RAW_R2"
    fi

    ## align read pairs
    BAM_FILE="${BAM_DIR}/${SAMPLE}.sorted.bam"
    if [[ -f "$BAM_FILE" ]]; then
	echo "[·] BAM exists for $SAMPLE."
    else
	echo "[*] Aligning $SAMPLE with BWAmeth"
	bwameth.py --reference "$REF_FASTA" --threads $THREADS \
            "$TRIMMED_R1" "$TRIMMED_R2" | \
	    samtools sort -@ $THREADS -o "$BAM_FILE" - && \
	    samtools index "$BAM_FILE"
    fi
done

echo '======================================'
echo 'All samples trimmed, aligned, indexed!'

