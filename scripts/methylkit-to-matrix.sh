#!/usr/bin/env bash

## generate methylation matrix from methylkit files
Rscript methylation-matrix.r

## restrict to illumina methylation array data
Rscript illumina-matrix.r

## estimate cell counts from methylation data
Rscript estimate-cell-counts.r

## generate DNA methylation scores
Rscript dna-methylation-scores.r

