#!/usr/bin/R

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
meth_file <- args[1]
output_file <- args[2]

meth_df <- data.frame(fread(meth_file)) 
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k$loc <- paste(ann450k$chr, ann450k$pos)    
meth_df$loc <- paste(meth_df$chr, meth_df$end)
merged_data <- merge(meth_df, ann450k, by = c("loc"))
colnames(merged_data)[colnames(merged_data) == "chr.x"] <- "chr"
methylation <- merged_data[, c("Name", grep("^(chr|start|end|X)", names(merged_data), value = TRUE))]
Illumina_matrix <- methylation[, names(methylation) != "chr.y"]
colnames(Illumina_matrix)[1] <- "CpGs"
write.csv(Illumina_matrix, file = output_file)
