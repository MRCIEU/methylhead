#!/usr/bin/env Rscript

library(methylKit)

args <- commandArgs(trailingOnly = TRUE)

samples <- args[1]
output_file <- args[2]
output_file2 <- args[3]

file <-read.csv(samples,header=F)
file<-data.frame(file)
file_paths <- file$V2
cleaned_paths <- gsub("\\[|\\]", "", file_paths)
cleaned_paths <- trimws(cleaned_paths)
    file.vector <- cleaned_paths
    file.vector <- file.vector[!grepl("control", file.vector)]
    file.vector <- file.vector[sapply(file.vector, function(f) 
       {
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
coverage_matrix <- meth_df[, grep("^(chr|start|end|coverage)", colnames(meth_df), value = TRUE)]
sample.ids <- grep("^(coverage|numCs|numTs|strand)", colnames(meth_df), invert = TRUE, value = TRUE)
coverage_matrix <-setNames(coverage_matrix, sample.ids)
columns_to_remove <- grep("^(coverage|numCs|numTs|strand)", colnames(meth_df), value = TRUE)
metyhlation_matrix <- meth_df[, setdiff(colnames(meth_df), columns_to_remove)]

write.csv(coverage_matrix,file = output_file2,row.names = FALSE)
write.csv(metyhlation_matrix, file = output_file, row.names = FALSE)
