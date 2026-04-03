#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

samples_file =args[1]  ## list of sample files with entropy data (csv format)
regions_file=args[2]  ## list of regions with frequencies (csv format)
entropy_file =args[3]  ## output file for methylation entropy matrix (csv format)
meth_file    =args[4]  ## output file for methylation level matrix (csv format)
sites_file   =args[5]  ## output file for number sites measured in each feature (csv format)
cores        =as.integer(args[6])

mincov=10
minsamples=0.5

library(parallel)
options(mc.cores=cores)

library(data.table)

samples = fread(samples_file)
regions = fread(regions_file)
colnames(regions) = c("chr", "start", "end",
regions = regions[regions$count > nrow(samples)*minsamples,]
regions$loc = paste(regions$chr, regions$start, sep=":")

ret = mclapply(samples$filename, function(filename) {
  cat(date(), " loading data from ", basename(filename),
      " (", which(filename == samples$filename), " of ", nrow(samples),
      ")\n")
  dat = fread(filename)
  dat$loc = paste(dat$chr, dat$start, sep=":")
  idx = match(regions$loc, dat$loc)
  dat = dat[idx,]
  c(dat$flip_pct, dat$meth_pct, dat$nmeth+dat$nunmeth)
})
ret = do.call(cbind, ret)
entropy = ret[1:nrow(regions),]
meth = ret[nrow(regions) + 1:nrow(regions),]
sites = ret[nrow(regions)*2 + 1:nrow(regions),]
colnames(entropy) = colnames(meth) = colnames(sites) = samples$sample_id

regions = regions[,c("chr","start","end")]

fwrite(cbind(regions, entropy), file=entropy_file, row.names = FALSE)
fwrite(cbind(regions, meth), file=meth_file, row.names = FALSE)
fwrite(cbind(regions, sites), file=sites_file, row.names = FALSE)
