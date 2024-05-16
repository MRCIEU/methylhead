#!/usr/bin/env Rscript
library(methylKit)
library(dplyr)
library(data.table)

file.vector <- list.files(pattern = "\\.markdup_CpG\\.methylKit", full.names = TRUE)
file.vector <- file.vector[!grepl("ctrl", file.vector)]
sample.ids <- gsub(".markdup_CpG.methylKit", "", basename(file.vector))
sample.ids<-as.list(sample.ids)
file.list<-as.list(file.vector)

myobj <- methRead(file.list,
           sample.id=sample.ids,
           pipeline = "amp",                
           assembly="hg19",
           treatment=c(rep(0,length(sample.ids))),
           mincov = 10)
myobj.filt <- filterByCoverage(myobj,
                      lo.count=10,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)

meth <- unite(myobj.filt, destrand=FALSE)
pm=percMethylation(meth)
pm=pm/100
meth_df<-data.frame(meth)
meth_df<-cbind(meth_df,(pm))
args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[1]
write.csv(meth_df,file = paste0(output_dir,"/DNAm_Full_Matrix.csv"),row.names=FALSE)
