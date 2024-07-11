#!/usr/bin/env Rscript

library(dplyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
pipeline <- args[1]
data_dir <- args[2]
dnam_dir <- args[3]
output_dir <- args[4]

if (pipeline == "bismark") {
    
    meth_df <- fread(paste0(dnam_dir,"/DNAm_Full_Matrix.csv"))
    meth_df <- data.frame(meth_df)
    
   
    blood_cell_types <- fread(paste0(data_dir,"/blood_cell_types.csv"))
    blood_cell_types <- data.frame(blood_cell_types)
    
    names(blood_cell_types) <- gsub("GSM\\d+_?", "", names(blood_cell_types))
    names(blood_cell_types) <- gsub(".Z\\d+_?", "", names(blood_cell_types))
    names(meth_df) <- gsub("^X", "", names(meth_df))
    names(meth_df) <- gsub("_R1_001_val_1\\d*", "", names(meth_df))
    
    sample.ids <- grep("_L001", names(meth_df), value = TRUE)
    regions <- blood_cell_types[c("chr","start","end")]
    
    CD4T <- rowMeans(blood_cell_types[c("Blood.T.CD4TT","Blood.T.CD4U7","Blood.T.CD4UM")])
    CD8T <- rowMeans(blood_cell_types[c("Blood.T.CD8TR","Blood.T.CD8U5","Blood.T.CD8UK")])
    NK <- rowMeans(blood_cell_types[c("Blood.NKTM","Blood.NKU1","Blood.NKUF")])
    Mono <- rowMeans(blood_cell_types[c("Blood.MonocytesTP","Blood.MonocytesU3","Blood.MonocytesUH")])
    Bcell <- rowMeans(blood_cell_types[c("Blood.BTX","Blood.BUB","Blood.BUR")])
    Granulocytes <- rowMeans(blood_cell_types[c("Blood.GranulocytesTZ","Blood.GranulocytesUD","Blood.GranulocytesUT")])
    beta.cell.types <- cbind(CD4T, CD8T, NK, Mono, Bcell, Granulocytes)        
    
    beta <- matrix(ncol=length(sample.ids), nrow=nrow(regions))
    colnames(beta) <- sample.ids
    colnames(beta) <- gsub("_L001\\d*", "", colnames(beta))
    
    calculate_mean_values <- function(data_frame, chr_value, start_value, end_value) {
      selected_data <- data_frame %>%
        filter(chr == chr_value & start >= start_value & end <= end_value)  
      selected_data <- selected_data[, -c(1:3)]  
      mean_values <- colMeans(selected_data) 
      return(mean_values)
    }
    
    meth_df2 <- meth_df[c("chr", "start", "end", unlist(sample.ids))]
    for (i in 1:nrow(regions)) {
      mean_values <- calculate_mean_values(meth_df2, regions$chr[i], regions$start[i], regions$end[i])
      beta[i, ] <- mean_values
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
    
    estimate_cell_counts <- estimate.cell.counts(beta, beta.cell.types)
    estimate_cell_counts_Normalized <- estimate_cell_counts / colSums(estimate_cell_counts)
    write.csv(estimate_cell_counts_Normalized, file = paste0(output_dir,"/estimate_cell_counts.csv"))
} else if (pipeline == "picard") {
   
   meth_df <- fread(paste0(dnam_dir,"/DNAm_Full_Matrix.csv"))
   meth_df <- data.frame(meth_df)
    
    args <- commandArgs(trailingOnly = TRUE)
    
    blood_cell_types <- fread(paste0(data_dir, "/blood_cell_types.csv"))
    blood_cell_types<-data.frame(blood_cell_types)

  names(blood_cell_types) <- gsub("GSM\\d+_?", "", names(blood_cell_types))
  names(blood_cell_types) <- gsub(".Z\\d+_?", "", names(blood_cell_types))
  names(meth_df) <- gsub("^X", "", names(meth_df))
  names(meth_df) <- gsub("_R1_001_val_1\\d*", "", names(meth_df))

sample.ids <- grep("_L001", names(meth_df), value = TRUE)
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
colnames(beta) <- gsub("_L001\\d*", "", colnames(beta))

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
estimate_cell_counts<-estimate.cell.counts(beta, beta.cell.types)
estimate_cell_counts_Normalized<-estimate_cell_counts/colSums(estimate_cell_counts)
write.csv(estimate_cell_counts_Normalized,file = paste0(output_dir,"/estimate_cell_counts.csv"))
}