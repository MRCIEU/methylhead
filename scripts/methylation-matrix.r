#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

samples_file <- args[1]  ## list of methylKit output files
scripts_dir <- args[2]
meth_file <- args[3]  ## methylation matrix (csv format)
cov_file  <- args[4]  ## read depth matrix (csv format)

library(data.table)

source(file.path(scripts_dir, "assemble-methylkits.r"))

samples <- read.csv(samples_file)

ret <- assemble.methylkits(
  samples,
  mincov=10,
  minsamples=0.5*nrow(samples))

fwrite(ret$meth, file=meth_file, row.names = FALSE)
fwrite(ret$coverage, file=cov_file, row.names = FALSE)


## previously used
## methylKit::methRead,filterByCoverage,unite
## replaced with assemble.methylseq()
## because methylKit handled gzip input
## by unzipping all inputs to /tmp (and never deleting them)
## and unnecessarily loading all data into memory


