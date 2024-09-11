#!/usr/bin/env Rscript

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
library(parallel)
options(mc.cores = 24)


args <- commandArgs(trailingOnly = TRUE)

meth_file <- args[1] 
output_file <- args[2] 

meth_df <- data.frame(fread(meth_file)) 
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k$loc <- paste(ann450k$chr, ann450k$pos)    
meth_df$loc <- paste(meth_df$chr, meth_df$end)
merged_data <- merge(meth_df, ann450k, by = c("loc"))
colnames(merged_data)[colnames(merged_data) == "chr.x"] <- "chr"
methylation <- merged_data[, c("Name", grep("^(X)", names(merged_data), value = TRUE))]
methylation <- methylation[, names(methylation) != "chr.y"]
colnames(methylation)[1] <- "CpGs"
twist_meth <- methylation
twist_meth <- data.frame(twist_meth)
rownames(twist_meth) <- twist_meth$CpGs
names(twist_meth) <- gsub("^X", "", names(twist_meth))
twist_meth <- twist_meth[!names(twist_meth) %in% "CpGs"]
meth_score <- data.matrix(twist_meth)
models <- meffonym.models()
errors <- c()
for (model in models) {
  tryCatch({
    score <- try(meffonym.score(meth_score, model)$score, silent = TRUE)
    if (class(score) == "try-error") {
      errors <- c(errors, model)
    }
  }, error = function(e) {
    errors <- c(errors, model)
  }, verbose = FALSE)
}
to_remove <- errors
models <- setdiff(models, to_remove)
scores_list <- list()
for (model in models) {
  score <- meffonym.score(meth_score, model)$score
  scores_list[[paste0("", model)]] <- score
}
meth_scores_all <- as.data.frame(scores_list)
rownames(meth_scores_all)<-colnames(meth_score)
DNA_Methylation_Scores <- data.frame(t(meth_scores_all))
names(DNA_Methylation_Scores) <- gsub("^X", "", names(DNA_Methylation_Scores))
write.csv(DNA_Methylation_Scores, file = output_file)
           
                              
