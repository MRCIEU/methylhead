#!/bin/bash 

DATA_DIR=/projects/MRC-IEU/research/projects/icep2/wp3/004/working/data

GENOME=hg19

CELL_TYPE_REGIONS=../data/target-regions.bed
CELL_TYPE_REFERENCE=../data/blood_cell_types_extended.zip..............
CELL_TYPE_REFERENCE=../data/blood-cell-type-reference.csv.gz
PANEL=../input/panel.csv

GENOME_REFERENCE=$DATA_DIR/genome-references
ORIGINAL_DATASET=$DATA_DIR/pipeline-test-dataset/raw
ORIGINAL_ALIGNMENT=$DATA_DIR/pipeline-test-dataset/aligned-files/bam

mkdir -p input
mkdir -p data
mkdir -p data/genome-reference
mkdir -p data/raw
TEST_PANEL=input/panel.csv
TEST_PANEL_BED=input/panel_$GENOME.bed
TEST_GENOME_REFERENCE=data/genome-reference
TEST_DATASET=data/raw
TEST_CELL_TYPE_REFERENCE=data/blood-cell-type-reference.csv

# download dataset https://www.ebi.ac.uk/ena/browser/view/PRJNA730913
bash scripts/download-dataset.sh $ORIGINAL_DATASET

## download and index genome reference hg19
bash ../scripts/create-reference.sh $GENOME $GENOME_REFERENCE

## align downloaded dataset to hg19
bash scripts/align-dataset.sh \
     $GENOME_REFERENCE \
     $ORIGINAL_DATASET \
     $ORIGINAL_ALIGNMENT

## select test panel of target regions
bash scripts/select-test-panel.sh \
     $ORIGINAL_ALIGNMENT \
     $PANEL \
     $CELL_TYPE_REGIONS \
     $TEST_PANEL \
     $TEST_PANEL_BED

## create test dataset blood cell type reference
Rscript scripts/create-test-cell-type-reference.r \
    $CELL_TYPE_REFERENCE \
    $TEST_PANEL \
    $TEST_CELL_TYPE_REFERENCE

## create test dataset genome reference to cover capture panel regions
bash scripts/create-test-genome.sh \
     $GENOME_REFERENCE/${GENOME}.fa \
     $TEST_PANEL_BED \
     $TEST_GENOME_REFERENCE/test.fa

bash ../scripts/create-reference.sh test $TEST_GENOME_REFERENCE

## create test dataset fastq files from aligned reads within capture panel regions
bash scripts/create-test-fastq-files.sh \
     $ORIGINAL_ALIGNMENT \
     $TEST_PANEL_BED \
     $TEST_DATASET

