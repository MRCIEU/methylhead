#!/usr/bin/R

library(methylKit)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(stringr)
library(GenomicRanges)
library(genomation)
library(meffonym)
library(dplyr)
library(ggplot2)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
pipeline <- args[1]
input_dir <- args[2]
output_dir <- args[3]
setwd(input_dir)

process_methylation_data <- function(file.vector, sample.ids, pipeline) {
  if (pipeline == "bismark") {
    
    file.list <- as.list(file.vector)
    myobj <- methRead(file.list,
                      sample.id = sample.ids,
                      pipeline = "bismarkCoverage",
                      assembly = "hg19",
                      treatment = rep(0, length(sample.ids)),
                      mincov = 10)
  } else if (pipeline == "picard") {

    file.list <- as.list(file.vector)
    myobj <- methRead(file.list,
                      sample.id = sample.ids,
                      pipeline = "amp",
                      assembly = "hg19",
                      treatment = rep(0, length(sample.ids)),
                      mincov = 10)
  } else {
    stop("Unknown pipeline")
  }
  
  myobj.filt <- filterByCoverage(myobj,
                                 lo.count = 10,
                                 lo.perc = NULL,
                                 hi.count = NULL,
                                 hi.perc = 99.9)
  
  meth <- unite(myobj.filt, destrand = FALSE)
  pm <- percMethylation(meth)
  pm <- pm / 100
  meth_df <- data.frame(meth)
  ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ann450k$loc <- paste(ann450k$chr, ann450k$pos)
    
  if (pipeline == "bismark") {
    meth_df$loc <- paste(meth_df$chr, meth_df$start)
  } else if (pipeline == "picard") {
    meth_df$loc <- paste(meth_df$chr, meth_df$end)
  }  
  meth_df <- cbind(meth_df, pm)
  merged_data <- merge(meth_df, ann450k, by = c("loc"))
  sample.ids2 <- unlist(sample.ids)
  mm <- merged_data[c("Name", sample.ids2)]
  methylation <- na.omit(mm)
  colnames(methylation)[1] <- "CpGs"
  return(list(methylation = methylation, meth = meth))
}

if (pipeline == "bismark") {
  file.vector <- list.files(pattern = "bismark\\.cov\\.gz", full.names = FALSE)
  file.vector <- file.vector[!grepl("ctrl", file.vector)]
  sample.ids <- gsub("_bismark_bt2_pe.deduplicated.bismark.cov.gz", "", basename(file.vector))
  sample.ids <- as.list(sample.ids)
} else if (pipeline == "picard") {
  file.vector <- list.files(pattern = "\\.markdup_CpG\\.methylKit", full.names = FALSE)
  file.vector <- file.vector[!grepl("ctrl", file.vector)]
  sample.ids <- gsub(".markdup_CpG.methylKit", "", basename(file.vector))
  sample.ids <- as.list(sample.ids)
} else {
  stop("Unknown pipeline")
}

result <- process_methylation_data(file.vector, sample.ids, pipeline)
methylation <- result$methylation
meth <- result$meth

if (pipeline == "bismark") {
  write.csv(methylation, file = paste0(output_dir, "/Methylation_matrix.csv"), row.names = FALSE)
} else if (pipeline == "picard") {
  write.csv(methylation, file = paste0(output_dir, "/Methylation_matrix.csv"), row.names = FALSE)
}

pdf(paste0(output_dir, "/correlation_plot.pdf"))
getCorrelation(meth, plot = TRUE)
dev.off()

pdf(paste0(output_dir, "/cluster_plot.pdf"))
clusterSamples(meth, dist = "correlation", method = "ward", plot = TRUE)
dev.off()