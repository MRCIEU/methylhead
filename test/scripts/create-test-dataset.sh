#!/bin/bash 

GENOME=http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
CELL_COUNT_REGIONS=../data/blood-cell-types-regions.bed
ORIGINAL_PANEL=../input/panel.csv
TEST_PANEL=input/panel.csv
TEST_REFERENCE=reference
TEST_DATASET=dataset

# download dataset https://www.ebi.ac.uk/ena/browser/view/PRJNA730913
bash scripts/download-dataset.sh original-dataset 

## download and index genome reference hg19, output in original-reference/
bash scripts/create-reference.sh $GENOME original-reference

## align downloaded dataset to hg19, output in original-alignment/
bash scripts/align-dataset.sh \
     original-reference \
     original-dataset \
     original-alignment

## select test panel of about 10K capture regions, save in $TEST_PANEL
bash scripts/select-test-panel.sh \
     original-alignment \
     $CELL_COUNT_REGIONS \
     $ORIGINAL_PANEL \
     $TEST_PANEL 

## create test dataset fastq files from aligned reads within capture panel regions
bash scripts/create-test-fastq-files.sh \
     original-alignment \
     $TEST_PANEL \
     $TEST_DATASET

## create test dataset genome reference to cover capture panel regions
bash scripts/create-test-reference.sh \
     original-reference \
     $TEST_PANEL \
     $TEST_REFERENCE
  
rm -rf original-reference original-alignment original-reference
