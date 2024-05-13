library(methylKit)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(stringr)
library(GenomicRanges)
library(genomation)
library(meffonym)
library(dplyr)
library(ggplot2)
library(ggcorrplot)
library(data.table)
library(ggpubr)

file.vector <- list.files(pattern = "_sorted\\.markdup_CpG\\.methylKit", full.names = TRUE)
file.vector <- file.vector[!grepl("ctrl", file.vector)]
sample.ids <- gsub("_sambamba_sorted.markdup_CpG.methylKit", "", basename(file.vector))
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

names(meth_df) <- gsub("^X", "", names(meth_df))
names(meth_df) <- gsub("_R1_001_val_1\\d*", "", names(meth_df))

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[2]

blood_cell_types <- fread(paste0(data_dir, "/blood_cell_types.csv"))
blood_cell_types<-data.frame(blood_cell_types)


names(blood_cell_types) <- gsub("GSM\\d+_?", "", names(blood_cell_types))
names(blood_cell_types) <- gsub(".Z\\d+_?", "", names(blood_cell_types))

regions<-blood_cell_types[c("chr","start","end")]

CD4T=rowMeans(blood_cell_types[c("Blood.T.CD4TT","Blood.T.CD4U7","Blood.T.CD4UM")])
CD8T=rowMeans(blood_cell_types[c("Blood.T.CD8TR","Blood.T.CD8U5","Blood.T.CD8UK")])
NK=rowMeans(blood_cell_types[c("Blood.NKTM","Blood.NKU1","Blood.NKUF")])
Mono=rowMeans(blood_cell_types[c("Blood.MonocytesTP","Blood.MonocytesU3","Blood.MonocytesUH")])
Bcell=rowMeans(blood_cell_types[c("Blood.BTX","Blood.BUB","Blood.BUR")])
Granulocytes=rowMeans(blood_cell_types[c("Blood.GranulocytesTZ","Blood.GranulocytesUD","Blood.GranulocytesUT")])
beta.cell.types<-cbind(CD4T,CD8T,NK,Mono,Bcell,Granulocytes)        

beta<-matrix(ncol=length(sample.ids),nrow=nrow(regions))
colnames(beta)<-sample.ids

calculate_mean_values <- function(data_frame, chr_value, start_value, end_value) {
  selected_data <- data_frame %>%
    filter(chr == chr_value & start >= start_value & end <= end_value)  
  selected_data <- selected_data[, -c(1:3)]  
  mean_values <- colMeans(selected_data) 
  return(mean_values)
}

meth_df2<-meth_df[c("chr","start","end",unlist(sample.ids))]
for (i in 1:nrow(regions)) {
  mean_values <- calculate_mean_values(meth_df2, regions$chr[i], regions$start[i], regions$end[i])
  beta[i,] <- mean_values
}

estimate.cell.counts <- function(beta, beta.cell.types) {
    stopifnot(nrow(beta) == nrow(beta.cell.types))
    apply(beta, 2, estimate.cell.counts0, beta.cell.types = beta.cell.types)
}

estimate.cell.counts0 <- function(beta, beta.cell.types) {
    stopifnot(length(beta) == nrow(beta.cell.types))
    require(quadprog)
    I <- diag(ncol(beta.cell.types))
    zero <- rep(0, ncol(beta.cell.types))
    loc.idx <- which(!is.na(beta))
    bcT.x.bc <- t(beta.cell.types[loc.idx,]) %*% beta.cell.types[loc.idx,]
    bcT.x.b <- t(beta.cell.types[loc.idx,]) %*% matrix(beta[loc.idx], nrow=length(loc.idx))
    counts <- solve.QP(bcT.x.bc, bcT.x.b, I, zero)$sol
    names(counts) <- colnames(beta.cell.types)
    counts
}

args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[1]
estimate_cell_counts<-estimate.cell.counts(beta, beta.cell.types)
estimate_cell_counts_Normalized<-estimate_cell_counts/colSums(estimate_cell_counts)
write.csv(estimate_cell_counts_Normalized,file = paste0(output_dir,"/estimate_cell_counts.csv"))
