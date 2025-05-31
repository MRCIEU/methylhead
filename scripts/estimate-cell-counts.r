#!/usr/bin/env Rscript

library(dplyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

blood_file        <- args[1] ## blood cell type reference data file
meth_file         <- args[2] ## DNA methylation file
output_file       <- args[3] ## cell count estimates output file

meth_df           <- data.frame(fread(meth_file))
blood_cell_types  <- fread(blood_file)
blood_cell_types  <- data.frame(blood_cell_types)

names(blood_cell_types) <- gsub("GSM\\d+_?",       "", names(blood_cell_types))
names(blood_cell_types) <- gsub(".Z\\d+_?",        "", names(blood_cell_types))
names(meth_df)          <- gsub("^X",              "", names(meth_df))
names(meth_df)          <- gsub("_R1_001_val_1\\d*", "", names(meth_df)) 

sample.ids <- names(meth_df %>% select(-chr, -start, -end))
regions    <- blood_cell_types[c("chr","start","end")]

CD4T        <- rowMeans(blood_cell_types[c("Blood.T.CD4TT",       "Blood.T.CD4U7",       "Blood.T.CD4UM")])
CD8T        <- rowMeans(blood_cell_types[c("Blood.T.CD8TR",       "Blood.T.CD8U5",       "Blood.T.CD8UK")])
NK          <- rowMeans(blood_cell_types[c("Blood.NKTM",          "Blood.NKU1",          "Blood.NKUF")])
Mono        <- rowMeans(blood_cell_types[c("Blood.MonocytesTP",   "Blood.MonocytesU3",   "Blood.MonocytesUH")])
Bcell       <- rowMeans(blood_cell_types[c("Blood.BTX",           "Blood.BUB",           "Blood.BUR")])
Granulocytes<- rowMeans(blood_cell_types[c("Blood.GranulocytesTZ","Blood.GranulocytesUD","Blood.GranulocytesUT")])
beta.cell.types <- cbind(CD4T, CD8T, NK, Mono, Bcell, Granulocytes)

setDT(meth_df)
meth_dt <- meth_df[, c("chr", "start", "end", sample.ids), with = FALSE]
setkey(meth_dt, chr, start, end)
setDT(regions)[, region_id := .I]
setkey(regions, chr, start, end)

hits <- foverlaps(meth_dt, regions, nomatch = 0L)
beta_dt <- hits[
  , lapply(.SD, mean),
    by = region_id,
    .SDcols = sample.ids
]
beta_dt <- merge(
  regions[, .(region_id)],
  beta_dt,
  by   = "region_id",
  all.x = TRUE
)[order(region_id)][, !"region_id"]

beta <- as.matrix(beta_dt)
colnames(beta) <- gsub("_L001\\d*", "", colnames(beta))

estimate.cell.counts <- function(beta, beta.cell.types) {
  stopifnot(nrow(beta) == nrow(beta.cell.types))
  require(quadprog)
  apply(beta, 2, function(b) {
    stopifnot(length(b) == nrow(beta.cell.types))
    loc <- which(!is.na(b) & rowSums(is.na(beta.cell.types)) == 0)
    bcT.bc <- t(beta.cell.types[loc, ]) %*% beta.cell.types[loc, ]
    bcT.b  <- t(beta.cell.types[loc, ]) %*% b[loc]
    diag(bcT.bc) <- diag(bcT.bc) + 1e-8
    solve.QP(bcT.bc, bcT.b, diag(ncol(beta.cell.types)), rep(0, ncol(beta.cell.types)))$sol -> res
    names(res) <- colnames(beta.cell.types)
    res
  })
}

estimate_cell_counts          <- estimate.cell.counts(beta, beta.cell.types)
estimate_cell_counts_norm     <- estimate_cell_counts / colSums(estimate_cell_counts)
estimate_cell_counts_norm[estimate_cell_counts_norm < 0] <- 1e-9
colnames(estimate_cell_counts_norm) <- gsub("\\.", "-", colnames(estimate_cell_counts_norm))
write.csv(estimate_cell_counts_norm, file = output_file)
