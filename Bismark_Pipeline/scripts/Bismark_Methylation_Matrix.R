#!/usr/bin/env Rscript
library(methylKit)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(stringr)
library(GenomicRanges)
library(genomation)
library(dplyr)
library(ggplot2)
library(ggcorrplot)
library(data.table)
library(ggpubr)

file.vector <- list.files(pattern = "bismark\\.cov\\.gz", full.names = FALSE)
file.vector <- file.vector[!grepl("ctrl", file.vector)]
sample.ids <- gsub("_bismark_bt2_pe.deduplicated.bismark.cov.gz", "", basename(file.vector))
sample.ids<-as.list(sample.ids)
file.list<-as.list(file.vector)

myobj <- methRead(file.list,
           sample.id=sample.ids,
           pipeline = "bismarkCoverage",                
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
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k$loc = paste(ann450k$chr, ann450k$pos)
meth_df$loc = paste(meth_df$chr,meth_df$start)
meth_df<-cbind(meth_df,(pm))
merged_data <- merge(meth_df, ann450k, by = c("loc"))
sample.ids2<-unlist(sample.ids)
mm <- merged_data[c("Name", sample.ids2)]
methylation<-na.omit(mm)
colnames(methylation)[1]<-"CpGs"

args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[1]

write.csv(methylation,file = paste0(output_dir,"/Bismark_methylation_matrix.csv"),row.names=FALSE)

pdf(paste0(output_dir,"/correlation_plot.pdf"))
getCorrelation(meth,plot=TRUE)
dev.off()

pdf(paste0(output_dir,"/culster_plot.pdf"))
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
dev.off()


