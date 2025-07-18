#!/usr/bin/env bash

REPO_DIR="$(realpath $1)"
DATA_DIR="$(realpath $2)"
PANEL_BED="$(realpath $3)"
OUT_DIR="$(realpath $4)"

SCRIPTS_DIR=$REPO_DIR/scripts

if [ "$#" -gt 4 ]; then
    ASSEMBLY=$5
else
    ASSEMBLY=hg19
fi

if [ "$ASSEMBLY" != "hg19" ] && [ "$ASSEMBLY" != "hg38" ]; then 
    echo "Assembly (arg 4) should be hg19 or hg38"
    exit 1
fi    

## create samples spreadsheet
ls $DATA_DIR/*.methylKit.gz | sed 's/\(.*\)/,\1/' > $OUT_DIR/samples.csv

## generate methylation matrix from methylkit files
Rscript $SCRIPTS_DIR/methylation-matrix.r \
    $OUT_DIR/samples.csv \
    $OUT_DIR/methylation.csv \
    $OUT_DIR/coverage.csv

## restrict to illumina methylation array data
Rscript $SCRIPTS_DIR/illumina-matrix.r \
    $OUT_DIR/methylation.csv \
    $OUT_DIR/illumina.csv \
    $ASSEMBLY

## estimate cell counts from methylation data
Rscript $SCRIPTS_DIR/estimate-cell-counts.r \
    $REPO_DIR/data/blood-cell-type-reference.csv.gz \
    $OUT_DIR/methylation.csv \
    $OUT_DIR/cell-counts.csv

## generate DNA methylation scores
Rscript $SCRIPTS_DIR/dna-methylation-scores.r \
    $OUT_DIR/methylation.csv \
    $OUT_DIR/methylation-scores.csv \
    $OUT_DIR/methylation-scores-sites.csv \
    $ASSEMBLY

## generate report
$SCRIPTS_DIR/bed2csv.sh $PANEL_BED $OUT_DIR/panel.csv
PWD=$(pwd)
cd $OUT_DIR
cp $SCRIPTS_DIR/methylkit-to-matrix.qmd .
quarto render methylkit-to-matrix.qmd \
    --to html \
    --output qc.html \
    -P panel=panel.csv \
    -P coverage_matrix=coverage.csv \
    -P cell_counts=cell-counts.csv \
    -P methylation_matrix=methylation.csv 



