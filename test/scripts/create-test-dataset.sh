#!/usr/bin/env bash
set -euo pipefail

#############################################################################
# Create a small test dataset 
#
# USAGE
#   bash scripts/create-test-dataset.sh <DATA_DIR> <GENOME_DIR>
#############################################################################

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <DATA_DIR> <GENOME_DIR>" >&2
  exit 1
fi

DATA_DIR=$(readlink -f "$1")
GENOME_DIR=$(readlink -f "$2")

REPO_BASE_DIR=$(readlink -f ..)

##################################################
## container config
##################################################

WGBS_SIF="methylhead-pipeline_wgbs_analysis.sif"
WGBS_ORAS="oras://docker.io/onuroztornaci/methylhead-pipeline:wgbs_analysis"

R_SIF="methylhead-pipeline_r.sif"
R_ORAS='oras://docker.io/onuroztornaci/methylhead-pipeline:qc_container'

BED_SIF="methylhead-pipeline_bedtools.sif"
BED_ORAS="docker://quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0"

rm -f "$WGBS_SIF"
apptainer pull "$WGBS_SIF" "$WGBS_ORAS"

rm -f "$R_SIF"
apptainer pull "$R_SIF" "$R_ORAS"

rm -f "$BED_SIF"
apptainer pull "$BED_SIF" "$BED_ORAS"

container_exec() {
    apptainer exec \
	--bind "$DATA_DIR":"$DATA_DIR" \
	--bind "$GENOME_DIR":"$GENOME_DIR" \
	--bind "$REPO_BASE_DIR":"$REPO_BASE_DIR" \
	"$@"
}

wgbs_container_exec() {
    container_exec "$WGBS_SIF" "$@"
}

r_container_exec() {
    container_exec "$R_SIF" "$@"
}

bed_container_exec() {
    container_exec "$BED_SIF" "$@"
}

##################################################
## convenience variables
##################################################

GENOME=hg19

CELL_TYPE_REGIONS=$REPO_BASE_DIR/data/target-regions.bed
CELL_TYPE_REFERENCE=$REPO_BASE_DIR/data/blood-cell-type-reference.csv.gz
PANEL=$REPO_BASE_DIR/input/panel.csv

RAW_DIR=$DATA_DIR/raw
ALIGN_DIR=$DATA_DIR/aligned-files/bam

mkdir -p data
mkdir -p data/genome-reference
mkdir -p data/raw

SCRIPTS=$(readlink -f scripts)
TEST_DIR=$(readlink -f .)

##################################################
## download dataset https://www.ebi.ac.uk/ena/browser/view/PRJNA730913
##################################################
bash $SCRIPTS/download-dataset.sh $RAW_DIR

##################################################
## download and index genome reference hg19
##################################################
wgbs_container_exec \
    bash $SCRIPTS/create-reference.sh $GENOME $GENOME_DIR

##################################################
## align downloaded dataset to hg19
##################################################
wgbs_container_exec \
    bash $SCRIPTS/align-dataset.sh $GENOME_DIR $RAW_DIR $ALIGN_DIR

##################################################
## select test panel of target regions
##################################################
wgbs_container_exec \
    bash $SCRIPTS/count-reads.sh \
        $ALIGN_DIR \
        $REPO_BASE_DIR/data/target-regions.bed \
        $TEST_DIR/input/read-counts-$GENOME.bed
r_container_exec \
    Rscript $SCRIPTS/select-test-panel.r \
        $REPO_BASE_DIR/input/panel.csv \
        $TEST_DIR/input/read-counts-$GENOME.bed \
        $TEST_DIR/input/panel.csv \
        $TEST_DIR/input/panel-$GENOME.bed

##################################################
## create test dataset blood cell type reference
##################################################
r_container_exec \
    Rscript $SCRIPTS/create-test-cell-type-reference.r \
        $REPO_BASE_DIR/data/blood-cell-type-reference.csv.gz \
        $TEST_DIR/input/panel.csv \
        $TEST_DIR/data/blood-cell-type-reference.csv

##################################################
## create test genome reference 
##################################################
bed_container_exec \
    bash $SCRIPTS/create-test-genome.sh \
       $GENOME_DIR/${GENOME}.fa \
       $TEST_DIR/input/panel-$GENOME.bed \
       $TEST_DIR/data/genome-reference/test.fa

##################################################
## index test genome
##################################################
wgbs_container_exec \
    bash $SCRIPTS/create-reference.sh \
        $TEST_DIR/data/genome-reference/test.fa \
        $TEST_DIR/data/genome-reference/

##################################################
## create test dataset fastq files with reads that overlap test panel
##################################################
wgbs_container_exec \
    bash $SCRIPTS/create-test-fastq-files.sh \
        $ALIGN_DIR \
        $TEST_DIR/input/panel-$GENOME.bed \
        $TEST_DIR/data/raw

