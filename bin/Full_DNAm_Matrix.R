#!/usr/bin/env Rscript

library(methylKit)

args <- commandArgs(trailingOnly = TRUE)
pipeline <- args[1]
samples <- args[2]
output_file <- args[3]

file <-read.csv(samples,header=F)
file<-data.frame(file)
file_paths <- file$V2
cleaned_paths <- gsub("\\[|\\]", "", file_paths)
cleaned_paths <- trimws(cleaned_paths)


if (pipeline == "bismark") {
    file.vector <- cleaned_paths
    file.vector <- file.vector[!grepl("control", file.vector)]
    sample.ids <- gsub("_bismark_bt2_pe.deduplicated.bismark.cov.gz", "", basename(file.vector))
    file.vector <- file.vector[sapply(file.vector, function(f) {
  file.size <- file.info(f)$size
  file.size >= 100 * 1024
  })] 
    sample.ids <- as.list(sample.ids)
    file.list <- as.list(file.vector)
       
    myobj <- methRead(file.list,
                      sample.id = sample.ids,
                      pipeline = "bismarkCoverage",                
                      assembly = "hg19",
                      treatment = c(rep(0, length(sample.ids))),
                      mincov = 10)
} else if (pipeline == "picard") {
    file.vector <- cleaned_paths
    file.vector <- file.vector[!grepl("control", file.vector)]
    file.vector <- file.vector[sapply(file.vector, function(f) {
  file.size <- file.info(f)$size
  file.size >= 100 * 1024
  })]
    sample.ids <- gsub(".markdup_CpG.methylKit", "", basename(file.vector))
    sample.ids <- as.list(sample.ids)
    file.list <- as.list(file.vector)
    myobj <- methRead(file.list,
                      sample.id = sample.ids,
                      pipeline = "amp",                
                      assembly = "hg19",
                      treatment = c(rep(0, length(sample.ids))),
                      mincov = 10)
} else {dos
    stop("Unknown pipeline")
}
myobj.filt <- filterByCoverage(myobj,
                               lo.count = 10,
                               lo.perc = NULL,
                               hi.count = NULL,
                               hi.perc = 99.9)
meth <- unite(myobj.filt, destrand = TRUE,min.per.group=as.integer(length(sample.ids)*0.50))
pm <- percMethylation(meth)
pm <- pm / 100
meth_df <- data.frame(meth)
meth_df <- cbind(meth_df, pm)
DNAm_Full_Matrix <- meth_df
write.csv(DNAm_Full_Matrix, file = output_file, row.names = FALSE)
