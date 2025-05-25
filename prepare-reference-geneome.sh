#!/usr/bin/env bash

set -euo pipefail

# -- Settings ------------------------------------------------------------
REF_DIR=${1:-reference}
FA=hg19.fa
URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/${FA}.gz"

# -- Container settings--
SIF_NAME="methylhead-pipeline_wgbs_analysis.sif"
SIF_ORAS="oras://docker.io/onuroztornaci/methylhead-pipeline:wgbs_analysis"

# -- Workflow ------------------------------------------------------------
mkdir -p "$REF_DIR" && cd "$REF_DIR"

rm -f "$SIF_NAME"
apptainer pull "$SIF_NAME" "$SIF_ORAS"

# 1) Download and decompress hg19 (if absent)
[ -s "$FA" ] || { wget -qO- "$URL" | gunzip -c > "$FA"; }

# 2) Create FASTA index (.fai) with samtools (if absent)
[ -s "$FA.fai" ] || apptainer exec "$SIF_NAME" samtools faidx "$FA"

# 3) Create sequence dictionary (.dict) with Picard (if absent)
DICT="${FA%.fa}.dict"
[ -s "$DICT" ] || apptainer exec "$SIF_NAME" picard CreateSequenceDictionary \
    REFERENCE="$FA" OUTPUT="$DICT"

echo "Reference genome prepared in $(pwd)"
