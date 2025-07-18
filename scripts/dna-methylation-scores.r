#!/usr/bin/env Rscript

args         <- commandArgs(trailingOnly = TRUE)
meth_file    <- args[1]
output_file  <- args[2]
sites_file   <- args[3]

assembly <- ifelse(length(args) > 3, args[4], "hg19")

stopifnot(assembly %in% c("hg19","hg38"))

if (assembly == "hg19") {
    annotations = c(
        "IlluminaHumanMethylation450kanno.ilmn12.hg19",
        "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
} else {
    annotations = c(
        "IlluminaHumanMethylationEPICv2anno.20a1.hg38")
}

for (annotation in annotations)
    library(annotation, character.only=T)
library(data.table)
library(meffonym)

annot <- lapply(annotations, function(annotation) getAnnotation(get(annotation)))
cols <- colnames(annot[[1]])
for (i in 1:length(annot)) cols <- intersect(cols, colnames(annot[[i]]))
annot <- lapply(annot, function(annot) annot[,cols])
annot <- do.call(rbind, annot)
annot$Name <- sub("_.*", "", annot$Name)
sites <- unique(annot$Name)
annot <- annot[match(sites, annot$Name),]

meth_df <- data.frame(fread(meth_file)) 

annot$loc <- paste(annot$chr, annot$pos)    
meth_df$loc <- paste(meth_df$chr, meth_df$end)
merged_data <- merge(meth_df, annot, by = c("loc"))
colnames(merged_data)[colnames(merged_data) == "chr.x"] <- "chr"
exclude_cols <- c("chr", "start", "end", "loc")
sample.ids <- setdiff(names(meth_df), exclude_cols)
methylation <- merged_data[, c("Name", sample.ids)]
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
colnames(DNA_Methylation_Scores) <- gsub("\\.", "-", colnames(DNA_Methylation_Scores))
write.csv(DNA_Methylation_Scores, file = output_file)                 

site_stats <- as.data.frame(t(sapply(models, function(model) {
    ret <- meffonym.get.model(model)
    num_sites_used <- length(intersect(ret$vars, rownames(meth_score)))
    num_sites_model <- length(ret$vars)
    c(sites=num_sites_used, model=num_sites_model, pct=round(num_sites_used/num_sites_model*100))
})))
colnames(site_stats) <- gsub("\\.", "-", colnames(site_stats))
write.csv(site_stats, file=sites_file)
        
