#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# CONFIG
###############################################################################
APPTAINER_BIN=apptainer

PIPELINE_REMOTE="oras://docker.io/onuroztornaci/methylhead-pipeline:wgbs_analysis"
PIPE_SIF=$(mktemp -u --suffix .sif)

echo "[*] Pulling WGBS container:   $PIPELINE_REMOTE"
"$APPTAINER_BIN" pull "$PIPE_SIF" "$PIPELINE_REMOTE"

###############################################################################
# ARGUMENTS
###############################################################################
if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <REF_FASTA> <RAW_FASTQ_DIR> <OUT_DIR> [THREADS]" >&2
  rm -f "$PIPE_SIF"; exit 1
fi

REF_FASTA=$(readlink -f "$1")
RAW_DIR=$(readlink -f "$2")
OUT_DIR=$(readlink -f "$3")
THREADS="${4:-8}"

TRIM_DIR="${OUT_DIR}/trimmed"
BAM_DIR="${OUT_DIR}/bam"
mkdir -p "$TRIM_DIR" "$BAM_DIR"

###############################################################################
# Helper wrapper
###############################################################################
pipe_exec() {
  "$APPTAINER_BIN" exec \
    --bind "$REF_FASTA":"$REF_FASTA" \
    --bind "$RAW_DIR":"$RAW_DIR" \
    --bind "$TRIM_DIR":"$TRIM_DIR" \
    --bind "$BAM_DIR":"$BAM_DIR" \
    "$PIPE_SIF" "$@"
}

###############################################################################
# 1) BWAmeth index (one-off)
###############################################################################
if ! ls "${REF_FASTA}".*bwameth.c2t* &>/dev/null; then
  echo "[*] Building BWAmeth index"
  pipe_exec bwameth.py index "$REF_FASTA"
else
  echo "[·] BWAmeth index already present."
fi

###############################################################################
# 2) FASTQ ➜ Trim ➜ BWAmeth ➜ BAM
###############################################################################
mapfile -t raw_r1_list < <(find "$RAW_DIR" -type f \( \
      -name '*_R1*.fastq.gz' -o -name '*_R1*.fq.gz' \
      -o -name '*_1*.fastq.gz'  -o -name '*_1*.fq.gz' \
    \) | sort)

[[ ${#raw_r1_list[@]} -eq 0 ]] && { echo "(!) No FASTQs found"; exit 1; }

for r1 in "${raw_r1_list[@]}"; do
  if [[ "$r1" =~ _R1 ]]; then
    r2="${r1/_R1/_R2}"; sample=$(basename "$r1" | sed -E 's/_R1[^.]*\.f(ast)?q\.gz//')
  else
    r2="${r1/_1/_2}";   sample=$(basename "$r1" | sed -E 's/_1[^.]*\.f(ast)?q\.gz//')
  fi
  [[ -f "$r2" ]] || { echo "(!) R2 missing for $sample, skipping."; continue; }

  trimmed_r1="${TRIM_DIR}/${sample}_val_1.fq.gz"
  trimmed_r2="${TRIM_DIR}/${sample}_val_2.fq.gz"

  if [[ -s "$trimmed_r1" && -s "$trimmed_r2" ]]; then
    echo "[·] $sample already trimmed."
  else
    echo "[*] Trimming $sample"
    pipe_exec trim_galore --paired --gzip --cores "$THREADS" \
        --basename "$sample" --output_dir "$TRIM_DIR" "$r1" "$r2"
  fi

  bam_out="${BAM_DIR}/${sample}.sorted.bam"
  if [[ -f "$bam_out" ]]; then
    echo "[·] BAM exists for $sample."
  else
    echo "[*] Aligning $sample with BWAmeth"
    pipe_exec bash -c "
      set -eo pipefail
      bwameth.py --reference '$REF_FASTA' --threads $THREADS \
        '$trimmed_r1' '$trimmed_r2' | \
      samtools sort -@ $THREADS -o '$bam_out' - && \
      samtools index '$bam_out'
    "
  fi
done

echo '======================================'
echo 'All samples trimmed, aligned, indexed!'

###############################################################################
# Cleanup
###############################################################################
echo "[*] Removing temporary container file"
rm -f "$PIPE_SIF"
