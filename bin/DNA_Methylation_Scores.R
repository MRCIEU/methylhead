#!/usr/bin/R 

library(meffonym)
library(dplyr)
library(data.table)
library(parallel)
options(mc.cores = 24)


args <- commandArgs(trailingOnly = TRUE)
pipeline <- args[1]
input_dir <- args[2]
output_dir <- args[3]

if (pipeline == "bismark") {
twist_meth<-fread(paste0(input_dir ,"/Methylation_matrix.csv"))
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
write.csv(bismark_scores,file = paste0(output_dir,"/Bismark_scores.csv"))   
write.csv(smoking_scores,file = paste0(output_dir,"/Bismark_smoking_scores.csv"))
           
} else if (pipeline == "picard") {

twist_meth<-fread(paste0(input_dir,"/Methylation_matrix.csv"))
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

write.csv(picard_scores,file = paste0(output_dir,"/Picard_scores.csv"))   
write.csv(smoking_scores,file = paste0(output_dir,"/Picard_smoking_scores.csv"))
}
           
                              
