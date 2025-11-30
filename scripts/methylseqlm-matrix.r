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

var.counts <- apply(counts,2,var,na.rm=T)
if (!all(is.na(var.counts))) {
  if (any(is.na(var.counts))) 
    cell.types <- colnames(counts)[!is.na(var.counts)]
  else
    cell.types <- colnames(counts)[-which.min(var.counts)]
} else
  cell.types <- c()

model.def <- readLines(model.file)

design = model.matrix(as.formula(model.def), pheno)
meth <- fread(meth.file, data.table=F)
coverage <- fread(coverage.file, data.table=F)

sites <- meth[,c("chr","start","end")]

common <- intersect(rownames(design), colnames(meth))

if (length(common) < ncol(meth))
  warning(
    ncol(meth)-length(common),
    " samples lost due to missing values in phenotype data")

meth.std = t(scale(t(meth[,common])))
meth.std = impute.fun(meth.std,1,mean,na.rm=T)
if (any(is.na(meth.std))) meth.std[is.na(meth.std)] = 0

segments <- methylseqlm.segment(
  meth.std,
  sites$chr,
  sites$start,
  design[common,], 
  max_dmr=20000,
  n0=1,m0=10,alpha=2)

common = setdiff(colnames(meth),colnames(sites))
meth.seqlm = methylseqlm.merge(
  as.matrix(meth[,common]),
  segments,
  as.matrix(coverage[,common]))

fwrite(
  data.frame(segments, meth.seqlm, check.names=F),
  file=output.file, row.names=FALSE)

