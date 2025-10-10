#!/usr/bin/env Rscript

library(methylKit)

args <- commandArgs(trailingOnly = TRUE)

samples   <- args[1]  ## list of methylKit output files
meth_file <- args[2]  ## methylation matrix (csv format)
cov_file  <- args[3]  ## read depth matrix (csv format)
assembly  <- args[4]  ## genome assembly (e.g., hg19, hg38)

files <- read.csv(samples,header=F)
paths <- files[[2]]
paths <- gsub("\\[|\\]", "", paths)
paths <- trimws(paths)

sample.ids <- gsub(".markdup_CpG.methylKit.gz", "", basename(paths))

mr.ret <- methRead(
  location = as.list(paths),
  sample.id = as.list(sample.ids),
  pipeline = "amp",                
  assembly = "assembly",
  treatment = rep(0, length(sample.ids)),
  mincov = 10)

mr.ret.filt <- filterByCoverage(
  mr.ret,
  lo.count = 10,
  lo.perc = NULL,
  hi.count = NULL,
  hi.perc = 99.9)

stats <- unite(mr.ret.filt, destrand = TRUE,min.per.group=as.integer(length(sample.ids)*0.50))
meth <- percMethylation(stats)
meth <- meth / 100
stats <- data.frame(stats,check.names=F)
coverage <- stats[,grep("coverage",colnames(stats))]
colnames(coverage) <- colnames(meth)
meth <- cbind(stats[,c("chr","start","end")], meth)
coverage <- cbind(stats[,c("chr","start","end")], coverage)

write.csv(coverage, file=meth_file, row.names = FALSE)
write.csv(meth, file=cov_file, row.names = FALSE)
