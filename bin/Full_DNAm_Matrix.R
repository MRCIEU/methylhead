#!/usr/bin/env Rscript

library(methylKit)
library(dplyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
pipeline <- args[1]
input_dir <- args[2]
output_dir <- args[3]


if (pipeline == "bismark") {
    file.vector <- list.files(path = input_dir, pattern = "bismark\\.cov\\.gz", full.names = FALSE)
    file.vector <- file.vector[!grepl("ctrl", file.vector)]
    sample.ids <- gsub("_bismark_bt2_pe.deduplicated.bismark.cov.gz", "", basename(file.vector))
    sample.ids <- as.list(sample.ids)
    file.list <- as.list(file.vector)

    myobj <- methRead(file.list,
                      sample.id = sample.ids,
                      pipeline = "bismarkCoverage",                
                      assembly = "hg19",
                      treatment = c(rep(0, length(sample.ids))),
                      mincov = 10)
} else if (pipeline == "picard") {
    file.vector <- list.files(path = input_dir, pattern = "\\.markdup_CpG\\.methylKit", full.names = TRUE)
    file.vector <- file.vector[!grepl("ctrl", file.vector)]
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

meth <- unite(myobj.filt, destrand = TRUE)
pm <- percMethylation(meth)
pm <- pm / 100
meth_df <- data.frame(meth)
meth_df <- cbind(meth_df, pm)

write.csv(meth_df, file = paste0(output_dir, "/DNAm_Full_Matrix.csv"), row.names = FALSE)
