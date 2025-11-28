library(data.table)

#' Merges methylkit output files for individual samples
#' into methylation and coverage matrices
#'
#' This is roughly like running methylKit::methRead()
#' methylKit::filterByCoverage() followed by methylKit::unite().
#' Unfortunately that sequence of calls is quite inflexible
#' for input files similar to but not exactly like methylkit, e.g.
#' CAMDA outputs (see camda-matrix.r). 
#' 
#' @param samples Data frame listing all samples identified by the 'sample_id'
#' column and corresponding methylkit output file identified by the 'filename'
#' column.
#' @param mincoverage Minimum coverage depth per sample
#' for each CpG site to be included
#' @param minsamples Minimum number of samples with minimum
#' coverage depth for each CpG site to be included
#' @param FUN Function for reading methylkit files 
#' @return List of two data frames, 'meth' and 'coverage'
#' each with the first three columns providing chromosome/start/end
#' and one column per sample. 'meth' provides DNAm levels
#' and 'coverage' provides number of reads on which DNAm levels
#' estimated.
assemble.methylkits <- function(
    samples,
    mincoverage,
    minsamples,
    FUN=function(filename) { data.table::fread(filename, data.table=F) }) {

  stopifnot(is.data.frame(samples))
  stopifnot(all(c("sample_id","filename") %in% colnames(samples)))
  
  required.columns <- c("chrBase","coverage","freqC","freqT")
  
  site_counts <- c()
  for (filename in samples$filename) {
    cat(date(), " loading sites from ", basename(filename),
        " (", which(filename == samples$filename), " of ", nrow(samples),
        ")\n")
    dat <- FUN(filename)
    stopifnot(all(required.columns %in% colnames(dat)))
    dat <- dat[dat$coverage >= mincoverage,]
    is.old <- names(site_counts) %in% dat$chrBase
    site_counts[is.old] <- site_counts[is.old] + 1 
    new.sites <- setdiff(dat$chrBase,names(site_counts))
    site_counts[new.sites] <- 1
  }		     
  
  site_counts = site_counts[site_counts > minsamples]
  
  sites <- data.frame(
    chrBase=names(site_counts),
    chr=sub("(chr[^\\.]+)\\.[0-9]+","\\1", names(site_counts)),
    base=as.integer(sub("chr[^\\.]+\\.([0-9]+)","\\1",names(site_counts))),
    count=unname(site_counts))
  
  sites <- sites[order(sites$chr,sites$base),]
  
  methcoverage <- sapply(samples$filename, function(filename) {
    cat(date(), " loading data from ", basename(filename),
        " (", which(filename == samples$filename), " of ", nrow(samples),
        ")\n")
    dat <- FUN(filename)	
    stopifnot(all(required.columns %in% colnames(dat)))
    dat <- dat[dat$coverage >= mincoverage,]
    dat$meth <- dat$freqC/(dat$freqC + dat$freqT)
    idx <- match(sites$chrBase, dat$chrBase)
    c(dat$meth[idx], dat$coverage[idx])
  })
  
  meth <- methcoverage[1:nrow(sites),]
  coverage <- methcoverage[nrow(sites) + 1:nrow(sites),]
  colnames(meth) <- colnames(coverage) <- samples$sample_id
  
  sites$start <- sites$base
  sites$end <- sites$base+1
  list(
    meth=cbind(sites[,c("chr","start","end")], meth),
    coverage=cbind(sites[,c("chr","start","end")], coverage))
}    
