#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

scripts.dir     <- args[1]
phenotype.file  <- args[2]
model.file    <- args[3]
meth.file       <- args[4]
counts.file     <- args[5]
coverage.file   <- args[6]
reads.file      <- args[7]
output.file      <- args[8]
cores <- as.integer(args[9])

library(data.table)
library(parallel)
library(impute)
options(mc.cores=cores)

source(file.path(scripts.dir, "methylseqlm.r"))
source(file.path(scripts.dir, "impute-function.r"))

## phenotypes ##
pheno <- data.frame(data.table::fread(phenotype.file))
stopifnot("sample_id" %in% colnames(pheno))

## read counts ##
reads <- data.table::fread(reads.file)
setnames(reads, "Reads passing filters", "reads")
reads$sample_id <- sub("_R[12]+$", "", reads$Sample)
reads$Sample <- NULL
reads = unique(reads)

pheno <- merge(pheno, reads, by = "sample_id")
rownames(pheno) <- pheno$sample_id

## cell counts ##
counts <- read.csv(counts.file,row.names=1,check.names=F)
counts <- t(counts)
pheno <- cbind(
  pheno,
  counts[match(pheno$sample_id,rownames(counts)),])

pheno$sample_id <- NULL

## remove cell counts that couldn't be estimated
var.counts <- apply(counts,2,var,na.rm=T)
if (!all(is.na(var.counts))) {
  if (any(is.na(var.counts))) 
    cell.types <- colnames(counts)[!is.na(var.counts)]
  else
    cell.types <- colnames(counts)[-which.min(var.counts)]
} else
  cell.types <- c()

## construct design matrix
model.def <- readLines(model.file)
old = options(na.action = "na.pass")
design = model.matrix(as.formula(model.def), pheno)
options(old)

## impute missing values in design matrix
if (any(is.infinite(design)))
  design[is.infinite(design)] = NA
if (any(is.nan(design)))
  design[is.nan(design)] = NA
is.bad.column = apply(design, 2, function(v) mean(is.na(v)) > 0.8)
is.bad.row = apply(design, 1, function(v) mean(is.na(v)) > 0.5)
if (any(is.bad.column))
  warning(sum(is.bad.column), " variables have too many missing values")
if (any(is.bad.row))
  warning(sum(is.bad.row), " samples have too many missing values")
if (nrow(design) < 10)
  warning(nrow(design), " is really too few samples for DMR analysis")
design = impute.knn(
  design[!is.bad.row,!is.bad.column],
  k=10, rowmax=0.5, colmax=0.8, maxp=1500, rng.seed=362436069
)$data

## remove variables with too little variance in design matrix
is.variable = apply(design,2,var,na.rm=T) > 0
pca = prcomp(design[,is.variable],center=T,scale=T)
is.variable = pca$sdev^2/sum(pca$sdev^2) > 0.01
design = cbind(intercept=1,pca$x[,is.variable])

## load methylation data
meth <- fread(meth.file, data.table=F)
coverage <- fread(coverage.file, data.table=F)
sites <- meth[,c("chr","start","end")]

## merge methylation and phenotype data
common <- intersect(rownames(design), colnames(meth))
if (length(common) < ncol(meth))
  warning(
    ncol(meth)-length(common),
    " methylation profiles lost due to missingness in phenotype data")

## scale methylation data
meth.std = t(scale(t(meth[,common])))
meth.std = impute.fun(meth.std,1,mean,na.rm=T)
if (any(is.na(meth.std))) meth.std[is.na(meth.std)] = 0

## segment methylome
segments <- methylseqlm.segment(
  meth.std, 
  sites$chr,
  sites$start,
  design[common,], 
  max_dmr=20000,
  n0=1,m0=10,alpha=2)

## calculate segment methylation levels
common = setdiff(colnames(meth),colnames(sites))
meth.seqlm = methylseqlm.merge(
  as.matrix(meth[,common]),
  segments,
  as.matrix(coverage[,common]))

## save outputs
fwrite(
  data.frame(segments, meth.seqlm, check.names=F),
  file=output.file, row.names=FALSE)

