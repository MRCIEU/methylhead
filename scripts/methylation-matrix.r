#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

samples_file=args[1]  ## list of methylKit output files (csv format)
sites_file  =args[2]  ## list of sites with frequencies (csv format)
meth_file   =args[3]  ## methylation matrix (csv format)
cov_file    =args[4]  ## read depth matrix (csv format)
cores       =as.integer(args[5])

mincov=10
minsamples=0.5

library(parallel)
options(mc.cores=cores)

library(data.table)

samples = fread(samples_file)
sites = fread(sites_file)
colnames(sites) = c("chr", "pos", "count")
sites = sites[sites$count > nrow(samples)*minsamples,]
sites$loc = paste(sites$chr, sites$pos, sep=":")

ret = mclapply(samples$filename, function(filename) {
  cat(date(), " loading data from ", basename(filename),
      " (", which(filename == samples$filename), " of ", nrow(samples),
      ")\n")
  dat = fread(filename)
  dat = dat[dat$coverage >= mincov,]
  dat$meth = dat$freqC/(dat$freqC + dat$freqT)  
  dat$loc = paste(dat$chr, dat$base, sep=":")
  idx = match(sites$loc, dat$loc)
  c(dat$meth[idx],dat$coverage[idx])
})
ret = do.call(cbind, ret)
meth = ret[1:nrow(sites),]
coverage = ret[nrow(sites) + 1:nrow(sites),]
colnames(meth) = colnames(coverage) = samples$sample_id

sites$start = sites$pos
sites$end = sites$pos+1
sites = sites[,c("chr","start","end")]

fwrite(cbind(sites, meth), file=meth_file, row.names = FALSE)
fwrite(cbind(sites, coverage), file=cov_file, row.names = FALSE)

