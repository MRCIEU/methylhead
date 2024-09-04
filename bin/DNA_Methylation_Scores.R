
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
library(parallel)
options(mc.cores = 24)
library(data.table)
library(parallel)
options(mc.cores = 24)


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
  
  meth <- unite(myobj.filt, destrand = TRUE)
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

if (pipeline == "bismark") {
twist_meth<-methylation
twist_meth<-data.frame(twist_meth)
rownames(twist_meth)<-twist_meth$CpGs
names(twist_meth) <- gsub("^X", "", names(twist_meth))
twist_meth<-twist_meth[!names(twist_meth)%in% "CpGs"]
bismark_score<-data.matrix(twist_meth)

models<-meffonym.models()
errors <- c()

for (model in models) {
  tryCatch({
    score <- try(meffonym.score(bismark_score, model)$score, silent = TRUE)
    if (class(score) == "try-error") {
      errors <- c(errors, model)
    }
  }, error = function(e) {
    errors <- c(errors, model)
  }, verbose = FALSE) 
}

to_remove<-errors
models <- setdiff(models, to_remove)

scores_list<-list()
for (model in models) {
  score <- meffonym.score(bismark_score, model)$score
  scores_list[[paste0("", model)]] <- score 
}   

bismark_scores_all<-as.data.frame(scores_list)
bismark_scores_all<-as.data.frame(scores_list)
bismark_scores<-data.frame(t(bismark_scores_all))
names(bismark_scores) <- gsub("^X", "", names(bismark_scores))
bismark_scores$Models<-models
smoking_scores <- bismark_scores[bismark_scores$Models %in% c("langdon-agnostic-current-vs-former", "langdon-agnostic-ever-vs-never", "langdon-candidate-current-vs-former","langdon-candidate-ever-vs-never"), ]

smoking_scores<-smoking_scores[!names(smoking_scores)%in%"Models"]
bismark_scores<-bismark_scores[!names(bismark_scores)%in%"Models"]
write.csv(bismark_scores,file = paste0(output_dir,"/DNA_Methylation_Scores.csv"))   
write.csv(smoking_scores,file = paste0(output_dir,"/Bismark_smoking_scores.csv"))
           
} else if (pipeline == "picard") {

twist_meth<- methylation
twist_meth<-data.frame(twist_meth)
rownames(twist_meth)<-twist_meth$CpGs
names(twist_meth) <- gsub("^X", "", names(twist_meth))
twist_meth<-twist_meth[!names(twist_meth)%in% "CpGs"]
picard_score<-data.matrix(twist_meth)

models<-meffonym.models()
errors <- c()

for (model in models) {
  tryCatch({
    score <- try(meffonym.score(picard_score, model)$score, silent = TRUE)
    if (class(score) == "try-error") {
      errors <- c(errors, model)
    }
  }, error = function(e) {
    errors <- c(errors, model)
  }, verbose = FALSE) 
}

to_remove<-errors
models <- setdiff(models, to_remove)

scores_list<-list()
for (model in models) {
  score <- meffonym.score(picard_score, model)$score
  scores_list[[paste0("", model)]] <- score 
}   

picard_scores_all<-as.data.frame(scores_list)
picard_scores_all<-as.data.frame(scores_list)
picard_scores<-data.frame(t(picard_scores_all))
names(picard_scores) <- gsub("^X", "", names(picard_scores))
picard_scores$Models<-models
smoking_scores <- picard_scores[picard_scores$Models %in% c("langdon-agnostic-current-vs-former", "langdon-agnostic-ever-vs-never", "langdon-candidate-current-vs-former","langdon-candidate-ever-vs-never"), ]

smoking_scores<-smoking_scores[!names(smoking_scores)%in%"Models"]
picard_scores<-picard_scores[!names(picard_scores)%in%"Models"]

write.csv(picard_scores,file = paste0(output_dir,"/DNA_Methylation_Scores.csv"))   
write.csv(smoking_scores,file = paste0(output_dir,"/Picard_smoking_scores.csv"))
}
           
                              
