#!/usr/bin/env Rscript

library(data.table)
library(quadprog)

args <- commandArgs(trailingOnly = TRUE)

ref_file <- args[1] ## blood cell type reference data file
meth_file <- args[2] ## DNA methylation file
output_file <- args[3] ## cell count estimates output file

meth_df <- fread(meth_file)
ref_df <- fread(ref_file, data.table=F)

coords <- c("chr","start","end")
sample.ids <- setdiff(colnames(meth_df), coords)
regions <- ref_df[,coords]

ref_df <- as.matrix(ref_df[,setdiff(colnames(ref_df),coords)])

setkey(meth_df, chr, start, end)
setDT(regions)[, region_id := .I]
setkey(regions, chr, start, end)

hits <- foverlaps(meth_df, regions, nomatch = 0L)

meth_df <- hits[
  , lapply(.SD, mean),
    by = region_id,
    .SDcols = sample.ids
]

meth_df <- merge(
  regions[, .(region_id)],
  meth_df,
  by = "region_id",
  all.x = TRUE
)[order(region_id)][, !"region_id"]

meth_df <- as.matrix(meth_df)

estimate.cell.counts <- function(meth, ref) {
  stopifnot(nrow(meth) == nrow(ref))
  apply(meth, 2, function(b) {
    stopifnot(length(b) == nrow(ref))
    loc <- which(!is.na(b) & rowSums(is.na(ref)) == 0)
    bcT.bc <- t(ref[loc, ]) %*% ref[loc, ]
    bcT.b <- t(ref[loc, ]) %*% b[loc]
    diag(bcT.bc) <- diag(bcT.bc) + 1e-8
    solve.QP(bcT.bc, bcT.b, diag(ncol(ref)), rep(0, ncol(ref)))$sol -> res
    names(res) <- colnames(ref)
    res
  })
}

counts <- estimate.cell.counts(meth_df, ref_df)
counts <- counts / colSums(counts)
counts[counts < 0] <- 1e-9

write.csv(counts, file = output_file)
