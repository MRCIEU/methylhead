#!/usr/bin/env Rscript

library(dplyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

cell_file      <- args[1] ## blood cell type reference data file
meth_file      <- args[2] ## DNA methylation file
output_file    <- args[3] ## cell count estimates output file

meth_df        <- data.frame(fread(meth_file))
cell_types     <- data.frame(fread(cell_file))

names(meth_df) <- gsub("^X",              "", names(meth_df))
names(meth_df) <- gsub("_R1_001_val_1\\d*", "", names(meth_df)) 

sample.ids     <- names(meth_df %>% select(-chr, -start, -end))
regions        <- cell_types[c("chr","start","end")]

combine_replicates <- function(df) {
  # Keep only data columns (exclude genomic coordinates)
  data_cols <- df[, setdiff(names(df), c("chr","start","end","startCpG","endCpG")), drop=FALSE]
  clean_names <- names(data_cols)
  # Remove GSM IDs (e.g. GSM1234567_SampleName → SampleName)
  clean_names <- gsub("^GSM[0-9]+[_\\.]", "", clean_names)
  # Normalize separators (replace '-' and '_' with '.', remove redundant dots)
  clean_names <- gsub("[-_]", ".", clean_names)
  clean_names <- gsub("\\.\\.", ".", clean_names)
  clean_names <- gsub("^\\.|\\.$", "", clean_names)
  # Blood T cell subtypes (map CenMem/EffMem/Naive/Eff → CD4 or CD8 from Lloyd et al. reference)
  clean_names <- gsub("CenMem[.-]?CD4", "CD4", clean_names)
  clean_names <- gsub("EffMem[.-]?CD4", "CD4", clean_names)
  clean_names <- gsub("Naive[.-]?CD4", "CD4", clean_names)
  clean_names <- gsub("Eff[.-]?CD4", "CD4", clean_names)
  clean_names <- gsub("CenMem[.-]?CD8", "CD8", clean_names)
  clean_names <- gsub("EffMem[.-]?CD8", "CD8", clean_names)
  clean_names <- gsub("Naive[.-]?CD8", "CD8", clean_names)
  clean_names <- gsub("Eff[.-]?CD8", "CD8", clean_names)
  # Split names into parts
  parts <- strsplit(clean_names, "\\.")
  # Whitelist of valid biological cell type tokens (from Lloyd et al. reference)
  whitelist <- c(
    "Blood","T","CD3","CD4","CD8","B","NK","Monocytes","Granulocytes",
    "Adipocytes","Aorta","Endothelium","Endocrine","Smooth","Muscle","Cardiomyocyte",
    "Fibroblasts","Macrophages","Osteoblasts","Podocytes",
    "Kidney","Glomerular","Tubular","Liver","Lung","Pancreas","Colon","Heart","Dermal","Skeletal",
    "Bladder","Prostate","Ovary","Endometrium","Fallopian","Thyroid","Tonsil","Tongue","Larynx","Esophagus","Pharynx","Breast",
    "Gastric","Small","int","Gallbladder","Pleura","Bone_marrow","Epidermal",
    "Neuron","Cortex","Cerebellum","Oligodendrocytes",
    "cfDNA","WBC"
  )
  # Remove non-informative suffixes (keep only meaningful tokens in whitelist)
  clean_parts <- lapply(parts, function(p) {
    while (length(p) > 1 && !(p[length(p)] %in% whitelist)) {
      p <- p[-length(p)]
    }
    p
  })
  # Recombine cleaned names
  clean_names <- sapply(clean_parts, function(p) paste(p, collapse="."))
  # Custom renaming rules
  clean_names <- gsub("^Blood\\.B$", "Bcell", clean_names)
  # Collapse replicates by group (rowMeans if multiple replicates)
  cols_by_group <- split(seq_along(clean_names), clean_names)
  collapsed <- sapply(cols_by_group, function(idx) {
    if (length(idx) == 1) data_cols[[idx]]
    else rowMeans(as.matrix(data_cols[, idx, drop=FALSE]))
  })
  return(as.matrix(collapsed))
}

beta.cell.types <- combine_replicates(cell_types)

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
