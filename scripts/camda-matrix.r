#!/usr/bin/env Rscript

library(data.table)
library(methylKit)

args <- commandArgs(trailingOnly = TRUE)
samples_file <- args[1] 
matrix_file <- args[2]
lower_file <- args[3]
upper_file <- args[4]

samples <- fread(samples_file, data.table=F)

outputs <- read.csv(text=sprintf("name,column,filename
matrix,ratio,%s
ci_upper,CI_lower,%s
ci_lower,CI_upper,%s",
matrix_file, lower_file, upper_file))
rownames(outputs) <- outputs$name

for (name in rownames(outputs))
  samples[[name]] <- ""

for (i in 1:nrow(samples)) {
  data <- data.table::fread(samples$filename[i])
  for (name in rownames(outputs)) {
    samples[[name]][i] <- file.path(
      dirname(samples$filename[i]),
      paste0(samples$sample_id[i], "_", name, ".tsv"))
    stat <- data[[outputs[name,"column"]]]
    fwrite(
      data.frame(
        chrBase=paste(data$chr,data$pos,sep="."),
        chr=data$chr,
        pos=data$pos,
        strand=ifelse(data$strand == "+", "F", "R"),
        coverage=11L,
        freqC=stat * 100,
        freqT=100-(stat * 100)),
      file=samples[[name]][i], sep="\t", quote=F)
  }
}

for (name in rownames(outputs)) {
  mret <- methylKit::methRead(
    location = as.list(samples[[name]]),
    sample.id = as.list(samples$sample_id),
    pipeline = "amp",                
    assembly = "assembly",
    treatment = rep(0,nrow(samples)),
    mincov = 10)
  
  mat <- methylKit::unite(
    mret,
    destrand = FALSE,
    min.per.group = as.integer(length(samples$sample_id) * 0.50))
  
  pm <-  methylKit::percMethylation(mat)
  pm <- pm / 100
  mat <- data.frame(mat)
  mat <- cbind(mat, pm)
  
  to_remove <- grep("(coverage|numCs|numTs|sample_id|strand)", colnames(mat))
  mat <- mat[,-to_remove]
  
  fwrite(mat, outputs[name,"filename"], row.names = FALSE)
  file.remove(samples[[name]])
}
