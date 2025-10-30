#!/usr/bin/env Rscript

library(data.table)
library(methylKit)

args <- commandArgs(trailingOnly = TRUE)

samples_file   <- args[1]  ## list of methylKit output files
meth_file <- args[2]  ## methylation matrix (csv format)
cov_file  <- args[3]  ## read depth matrix (csv format)

samples <- read.csv(samples_file)

mr.ret <- methylKit::methRead(
  location = as.list(samples$filename),
  sample.id = as.list(samples$sample_id),
  pipeline = "amp",                
  assembly = "assembly",
  treatment = rep(0, nrow(samples)),
  mincov = 10)

mr.ret.filt <- methylKit::filterByCoverage(
  mr.ret,
  lo.count = 10,
  lo.perc = NULL,
  hi.count = NULL,
  hi.perc = 99.9)

stats <- methylKit::unite(mr.ret.filt, destrand = TRUE,min.per.group=as.integer(nrow(samples)*0.50))
meth <- methylKit::percMethylation(stats)
meth <- meth / 100
stats <- data.frame(stats,check.names=F)
coverage <- stats[,grep("coverage",colnames(stats))]
colnames(coverage) <- colnames(meth)
meth <- cbind(stats[,c("chr","start","end")], meth)
coverage <- cbind(stats[,c("chr","start","end")], coverage)

fwrite(meth, file=meth_file, row.names = FALSE)
fwrite(coverage, file=cov_file, row.names = FALSE)
