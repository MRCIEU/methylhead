#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Download Docker image (if not already present)
# ------------------------------------------------------------

APPTAINER_BIN=apptainer
CONTAINER_REMOTE="oras://docker.io/onuroztornaci/methylhead-pipeline:wgbs_analysis"
TMP_CONTAINER=$(mktemp -u --suffix .sif)

# Pull the container image only once
echo "[*] Pulling Apptainer container: $CONTAINER_REMOTE"
"$APPTAINER_BIN" pull "$TMP_CONTAINER" "$CONTAINER_REMOTE"

#############################################################
#  End-to-end WGBS pipeline: Trim Galore!  ➜  BWAmeth  ➜  sorted+indexed BAM (with Docker)
#
#  Usage:
#    ./align-dataset.sh <REF_FASTA> <RAW_FASTQ_DIR> <OUT_DIR> [THREADS]
#
#  Outputs:
#    $OUT_DIR/trimmed/   → *_val_1.fq.gz / *_val_2.fq.gz
#    $OUT_DIR/bam/       → <sample>.sorted.bam + <sample>.sorted.bam.bai
#
#  All commands are run inside the $DOCKER_IMAGE container.
#############################################################

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <REF_FASTA> <RAW_FASTQ_DIR> <OUT_DIR> [THREADS]" >&2
  rm -f "$TMP_CONTAINER"
  exit 1
fi

REF_FASTA=$(readlink -f "$1")
RAW_DIR=$(readlink -f "$2")
OUT_DIR=$(readlink -f "$3")
THREADS="${4:-8}"

TRIM_DIR="${OUT_DIR}/trimmed"
BAM_DIR="${OUT_DIR}/bam"

mkdir -p "$TRIM_DIR" "$BAM_DIR"

# Find raw FASTQ pairs (R1 only)
mapfile -t raw_r1_list < <(find "$RAW_DIR" -type f \( \
      -name '*_R1*.fastq.gz' -o -name '*_R1*.fq.gz' \
      -o -name '*_1*.fastq.gz'  -o -name '*_1*.fq.gz' \
    \) | sort)

if [[ ${#raw_r1_list[@]} -eq 0 ]]; then
  echo "(!) No R1 fastq files found in ${RAW_DIR}" >&2
  rm -f "$TMP_CONTAINER"
  exit 1
fi

apptainer_exec() {
  # Wrapper for running commands inside the container with all folders bound
  # Usage: apptainer_exec CMD...
  "$APPTAINER_BIN" exec \
    --bind "$REF_FASTA":"$REF_FASTA" \
    --bind "$RAW_DIR":"$RAW_DIR" \
    --bind "$TRIM_DIR":"$TRIM_DIR" \
    --bind "$BAM_DIR":"$BAM_DIR" \
    "$TMP_CONTAINER" "$@"
}

for r1 in "${raw_r1_list[@]}"; do
  if [[ "$r1" =~ _R1 ]]; then
    r2="${r1/_R1/_R2}"
    sample=$(basename "$r1" | sed -E 's/_R1[^.]*\.f(ast)?q\.gz//')
  else
    r2="${r1/_1/_2}"
    sample=$(basename "$r1" | sed -E 's/_1[^.]*\.f(ast)?q\.gz//')
  fi

  if [[ ! -f "$r2" ]]; then
    echo "(!) Mate R2 not found for $r1, skipping." >&2
    continue
  fi

  trimmed_r1="${TRIM_DIR}/$(basename "$r1" | sed -E 's/\.f(ast)?q\.gz$/_val_1.fq.gz/')"
  trimmed_r2="${trimmed_r1/_val_1/_val_2}"

  # Trimming
  if [[ -s "$trimmed_r1" && -s "$trimmed_r2" ]]; then
    echo "[·] $sample already trimmed, skipping."
  else
    echo "[*] Trimming: $sample"
    apptainer_exec trim_galore --paired --gzip --cores "$THREADS" \
      --output_dir "$TRIM_DIR" \
      "$r1" "$r2"
  fi

  # Alignment
  if [[ ! -s "$trimmed_r1" || ! -s "$trimmed_r2" ]]; then
    echo "(!) Trimmed fastq missing for $sample, alignment skipped." >&2
    continue
  fi

  bam_out="${BAM_DIR}/${sample}.sorted.bam"
  if [[ -f "$bam_out" ]]; then
    echo "[·] BAM already exists for $sample, skipping alignment."
    continue
  fi

  echo "[*] BWAmeth alignment: $sample"
  apptainer_exec bash -c "
    bwameth.py --reference '$REF_FASTA' --threads $THREADS \
      '$trimmed_r1' '$trimmed_r2' | \
    samtools sort -@ $THREADS -o '$bam_out' - && \
    samtools index '$bam_out'
  "
  echo "[+] Done: $(basename "$bam_out")"
done

echo "========================================"
echo "All pairs trimmed, aligned, and indexed!"

# Cleanup container image
echo "[*] Cleaning up container file: $TMP_CONTAINER"
rm -f "$TMP_CONTAINER"