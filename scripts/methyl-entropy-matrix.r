#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

samples_file=args[1]  ## sample files with entropy data (csv format)
regions_file=args[2]  ## list of regions with frequencies (csv format)
flip_file   =args[3]  ## output methylation flips matrix (csv format)
entropy_file=args[4]  ## output methylation entropy matrix (csv format)
meth_file   =args[5]  ## output methylation level matrix (csv format)
sites_file  =args[6]  ## output number sites measured (csv format)
cores        =as.integer(args[7])

mincov=10
minsamples=0.5

library(parallel)
options(mc.cores=cores)

library(data.table)

samples = fread(samples_file)
regions = fread(regions_file)
colnames(regions) = c("chr", "start", "end", "count")
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

  ## conditional entropy
  pu = dat$nmeth/(dat$nmeth + dat$nunmeth) ## P(unmeth)
  pm = dat$nunmeth/(dat$nmeth + dat$nunmeth) ## P(meth)
  puu = dat$nuu/(dat$nuu + dat$num) ## P(unmeth|unmeth before)
  pmu = dat$num/(dat$nuu + dat$num) ## P(meth|unmeth before)
  pum = dat$nmu/(dat$nmm + dat$nmu) ## P(unmeth|meth before)
  pmm = dat$nmm/(dat$nmm + dat$nmu) ## P(meth|meth before)
  dat$entropy = (
    -pu*(puu*log2(puu) + pmu*log2(pmu))
    -pm*(pum*log2(pum) + pmm*log2(pmm))) 

  c(dat$flip_pct, dat$entropy, dat$meth_pct, dat$nmeth+dat$nunmeth)
})
ret = do.call(cbind, ret)
flip = ret[1:nrow(regions),]
entropy = ret[nrow(regions) + 1:nrow(regions),]
meth = ret[2*nrow(regions) + 1:nrow(regions),]
sites = ret[3*nrow(regions) + 1:nrow(regions),]
colnames(flip) = colnames(entropy) = colnames(meth) = colnames(sites) = samples$sample_id

regions = regions[,c("chr","start","end")]

fwrite(cbind(regions, flip), file=flip_file, row.names = FALSE)
fwrite(cbind(regions, entropy), file=entropy_file, row.names = FALSE)
fwrite(cbind(regions, meth), file=meth_file, row.names = FALSE)
fwrite(cbind(regions, sites), file=sites_file, row.names = FALSE)
