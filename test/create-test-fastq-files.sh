#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Region-targeted BAM ➜ FASTQ ➜ BWAmeth ➜ CpG (Apptainer only)
#
# USAGE
#   ./create-test-fastq-files.sh <BAM_DIR> <TARGETS.bed> <OUT_DIR> <REF.fa> \
#                                [THREADS=8] [PAD=100]
#
# WHAT IT DOES
#   0. Pull two SIF images (wgbs_analysis + bedtools) once, then cache.
#   1. Create a padded BED (+PAD bp) once with Bedtools.
#   2. For every *.bam:
#        • filter reads to BED, export paired FASTQs
#        • BWAmeth realignment → coordinate-sorted BAM
#        • MethylDackel CpG extraction on targets
#   3. Cleans heavy BAMs to save space.
#   All tools run inside containers; no sudo / local installs needed.
###############################################################################

# -------- CONFIG -------- #
APPTAINER=${APPTAINER_BIN:-apptainer}

WGBS_REMOTE="oras://docker.io/onuroztornaci/methylhead-pipeline:wgbs_analysis"
BED_REMOTE="docker://quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0"
THREADS="${5:-8}"
PAD="${6:-100}"
# ------------------------ #

[[ $# -lt 4 ]] && {
  echo "Usage: $0 <BAM_DIR> <TARGETS.bed> <OUT_DIR> <REF.fa> [THREADS] [PAD]" >&2
  exit 1
}

BAM_DIR=$(readlink -f "$1")
ORIG_BED=$(readlink -f "$2")
OUT_DIR=$(readlink -f "$3")
REF_FA=$(readlink -f "$4")
REF_DIR=$(dirname "$REF_FA")        # bind the whole folder
REF_BAS=$(basename "$REF_FA")

mkdir -p "$OUT_DIR"

# -------- Pull images (cached in OUT_DIR) -------- #
WGBS_SIF="$OUT_DIR/wgbs_analysis.sif"
BED_SIF="$OUT_DIR/bedtools.sif"

[[ -s $WGBS_SIF ]] || {
  echo "[*] Pulling wgbs_analysis image…"
  $APPTAINER pull "$WGBS_SIF" "$WGBS_REMOTE"
}
[[ -s $BED_SIF ]] || {
  echo "[*] Pulling Bedtools image…"
  $APPTAINER pull "$BED_SIF" "$BED_REMOTE"
}

# -------- Build padded BED once -------- #
PAD_BED="$OUT_DIR/regions_padded_${PAD}bp.bed"
if [[ ! -s $PAD_BED ]]; then
  echo "[*] Building padded BED (±${PAD} bp)"
  $APPTAINER exec \
    --bind "$ORIG_BED":/in/targets.bed:ro \
    --bind "$REF_DIR":/ref:ro \
    "$BED_SIF" bash -euo pipefail -c "
      sort -k1,1 -k2,2n /in/targets.bed |
      bedtools merge -i - |
      bedtools slop -b $PAD -g /ref/${REF_BAS}.fai -i - > $PAD_BED.tmp
    "
  mv "$PAD_BED.tmp" "$PAD_BED"
fi

###############################################################################
# helper: run inside wgbs container with proper mounts
###############################################################################
wgbs() {
  $APPTAINER exec \
    --bind "$BAM_DIR":/bams:ro \
    --bind "$OUT_DIR":/out \
    --bind "$PAD_BED":/bed/targets.bed:ro \
    --bind "$REF_DIR":/ref:ro \
    "$WGBS_SIF" "$@"
}

###############################################################################
# main loop
###############################################################################
shopt -s nullglob
for BAM in "$BAM_DIR"/*.bam; do
  SAMPLE=$(basename "$BAM" | cut -d. -f1)
  echo -e "\n========= $SAMPLE ========="

  wgbs bash -euo pipefail -c "
    # 1. filter → name-sort
    samtools view -b -L /bed/targets.bed /bams/$(basename "$BAM") |
      samtools sort -n -@ $THREADS -o /out/${SAMPLE}.namesort.bam -

    # 2. FASTQ export
    samtools fastq -@ $THREADS \
      -1 /out/${SAMPLE}_1.fastq.gz \
      -2 /out/${SAMPLE}_2.fastq.gz \
      -0 /dev/null -s /dev/null /out/${SAMPLE}.namesort.bam
    rm /out/${SAMPLE}.namesort.bam

    # 3. BWAmeth realign
    bwameth.py --threads $THREADS --reference /ref/${REF_BAS} \
      /out/${SAMPLE}_1.fastq.gz /out/${SAMPLE}_2.fastq.gz |
      samtools sort -@ $THREADS -o /out/${SAMPLE}.aligned.bam -
    samtools index /out/${SAMPLE}.aligned.bam

    # 4. CpG calling (targets only)
    MethylDackel extract -@ $THREADS -l /bed/targets.bed \
      -o /out/${SAMPLE}_methylation \
      /ref/${REF_BAS} /out/${SAMPLE}.aligned.bam

    # 5. optional disk cleanup
    rm /out/${SAMPLE}.aligned.bam /out/${SAMPLE}.aligned.bam.bai
  "
done
shopt -u nullglob

do_cleanup() {
  echo -e "\n[*] Final cleanup…"

  # 1. Delete CpG bedGraphs (keep .txt summary if you like)
  find "$OUT_DIR" -type f -name '*_methylation_CpG.bedGraph' -print -delete

  # 2. Remove padded BED if you don’t care to reuse
  rm -f "$PAD_BED"

  # 3. Purge container images (comment out to keep cache)
  rm -f "$WGBS_SIF" "$BED_SIF"
}

# run even if script exits early
trap do_cleanup EXIT
echo -e "\n🎉  Pipeline finished — results in $OUT_DIR"
