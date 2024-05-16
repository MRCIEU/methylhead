#!/usr/bin/env Rscript
library(meffonym)
library(dplyr)
library(data.table)

twist_meth<-fread("picard_methylation.csv")
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

args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[1]

write.csv(picard_scores,file = paste0(output_dir,"/Picard_scores.csv"))   
write.csv(smoking_scores,file = paste0(output_dir,"/Picard_smoking_scores.csv"))
           
                              