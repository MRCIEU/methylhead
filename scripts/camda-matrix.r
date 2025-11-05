#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
samples_file <- args[1]
scripts_dir <- args[2]
matrix_file <- args[3]

library(data.table)

source(file.path(scripts_dir, "assemble-methylkits.r"))

samples <- fread(samples_file, data.table=F)

required.columns <- c("chr","pos","strand","ratio","CI_upper","CI_lower")

ret <- assemble.methylkits(
  samples,
  mincov=10,
  minsamples=0.5*nrow(samples),
    FUN=function(filename) {
      dat <- fread(filename, data.table=F)
      dat <- dat[which(dat$CI_upper - dat$CI_lower > 0.5),]
      stopifnot(all(required.columns %in% colnames(dat)))
      data.frame(
        chrBase=paste(dat$chr,dat$pos,sep="."),
        chr=dat$chr,
        base=dat$pos,
        strand=ifelse(dat$strand == "+", "F", "R"),
        coverage=11L,
        freqC=dat$ratio * 100,
        freqT=100-(dat$ratio * 100))
    })

fwrite(ret$meth, file=matrix_file, row.names = FALSE)

