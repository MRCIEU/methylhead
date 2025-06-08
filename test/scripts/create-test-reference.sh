#!/usr/bin/env bash
set -euo pipefail

################################################################################
# create-test-reference.sh   (container-aware)
#
# Usage:
#   bash create-test-reference.sh <reference_fasta> <regions_bed> \
#                                 <output_prefix> <output_dir> [EXT_BP]
#
#   <reference_fasta> : Whole-genome FASTA
#   <regions_bed>     : Original (un-padded) BED
#   <output_prefix>   : Prefix for all output files (e.g. test-ref)
#   <output_dir>      : Destination directory
#   [EXT_BP]          : Optional; bp to extend on both sides (default 0)
#   test-target-regions.bed: Filtered BED with intervals within chromosome bounds
#   which is subset of the original BED file and includes whole original regions.
################################################################################

###############################################################################
# 0. Container settings (leave the rest of the code unchanged)
###############################################################################
USE_CONTAINER=${USE_CONTAINER:-1}            # 1 = use SIF images, 0 = native

APPTAINER=${APPTAINER_BIN:-apptainer}
WGBS_IMG="oras://docker.io/onuroztornaci/methylhead-pipeline:wgbs_analysis"
BED_IMG="docker://quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CACHE_DIR="${SCRIPT_DIR}/.sif-cache"
mkdir -p "$CACHE_DIR"

WGBS_SIF="$CACHE_DIR/wgbs_analysis.sif"
BED_SIF="$CACHE_DIR/bedtools.sif"

if [[ $USE_CONTAINER -eq 1 ]]; then
  [[ -s $WGBS_SIF ]] || $APPTAINER pull "$WGBS_SIF" "$WGBS_IMG"
  [[ -s $BED_SIF  ]] || $APPTAINER pull "$BED_SIF"  "$BED_IMG"

  samtools()   { $APPTAINER exec --bind "$PWD":/work "$WGBS_SIF" samtools   "$@"; }
  bwameth.py() { $APPTAINER exec --bind "$PWD":/work "$WGBS_SIF" bwameth.py "$@"; }
  bedtools()   { $APPTAINER exec --bind "$PWD":/work "$BED_SIF"  bedtools   "$@"; }
fi
###############################################################################

if [[ $# -lt 4 || $# -gt 5 ]]; then
  echo "Usage: $0 <reference_fasta> <regions_bed> <output_prefix> <output_dir> [EXT_BP]"
  exit 1
fi

REF="$(readlink -f "$1")"
BED="$(readlink -f "$2")"
OUTPREFIX="$3"
OUTDIR="$(readlink -f "$4")"
EXT="${5:-0}"              # extension bp; default 0

mkdir -p "$OUTDIR"
cd "$OUTDIR"

###############################################################################
# 1. Extend BED ±EXT bp
###############################################################################
BED_EXT="${OUTPREFIX}_slop${EXT}.bed"

# Create REF.fai if missing (required by bedtools slop)
[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

bedtools slop -b "$EXT" -i "$BED" -g "${REF}.fai" \
  | awk '($2<0){$2=0}1' OFS='\t' > "$BED_EXT"

###############################################################################
# 2. Extract FASTA from the extended BED
###############################################################################
bedtools getfasta -fi "$REF" -bed "$BED_EXT" -fo regions_tmp.fa

###############################################################################
# 3. Merge all segments from the same chromosome into a single FASTA entry
###############################################################################
awk '
/^>/ {
  split($0, arr, "[:]");
  chr = substr(arr[1], 2);
  if (!seq[chr]) order[++n] = chr;
  curr_chr = chr;
  next;
}
{ seq[curr_chr] = seq[curr_chr] $0 }
END {
  for (i = 1; i <= n; i++) {
    print ">" order[i];
    s = seq[order[i]];
    for (j = 1; j <= length(s); j += 60)
      print substr(s, j, 60);
  }
}
' regions_tmp.fa > "${OUTPREFIX}.fa"
rm regions_tmp.fa

###############################################################################
# 4. FASTA index (.fai) and sequence dictionary (.dict)
###############################################################################
samtools faidx "${OUTPREFIX}.fa"  > /dev/null
samtools dict "${OUTPREFIX}.fa"   > "${OUTPREFIX}.dict"

###############################################################################
# 5. Filter the original BED so all intervals are within chromosome bounds
###############################################################################
awk 'NR==FNR {chrlen[$1]=$2; next} $3 <= chrlen[$1]' \
    "${OUTPREFIX}.fa.fai" "$BED" > test-target-regions.bed

###############################################################################
# 6. Build the BWAmeth index
###############################################################################
bwameth.py index "${OUTPREFIX}.fa"

###############################################################################
# 7. Summary
###############################################################################
do_cleanup() {
  echo -e "\n[*] Final cleanup…"
  # 1. Remove padded BED if you don’t care to reuse
  find "$OUT_DIR" -type -f -name "${OUTPREFIX}_slop${EXT}.bed" -print -delete
  # 2. Purge container images (comment out to keep cache)
  rm -f "$WGBS_SIF" "$BED_SIF"
}

echo -e "\nReference built and indexed in $OUTDIR:"
echo "- FASTA:            ${OUTPREFIX}.fa"
echo "- Index files:      ${OUTPREFIX}.fa.fai, ${OUTPREFIX}.dict, *bwameth.c2t.*"
echo "- Filtered BED:     test-target-regions.bed"

