#!/usr/bin/env R

### Usage ###
# Rscript test-dnam-assocs.r <data_dir> <phenotypes.csv> [<output_dir>]
# e.g.
#   Rscript test-dnam-assocs.r data/ phenotypes.csv output

library(data.table)
library(dplyr)
library(ewaff)

args <- commandArgs(trailingOnly = TRUE)

file_list_path  <- args[1]
phenotype.file  <- args[2]
models.file     <- args[3]
output.dir      <- args[4]

dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(output.dir, "sessionInfo.txt")); sessionInfo(); sink()

key_map <- c(
  "methylation-matrix.csv"            = "Methylation_matrix",
  "illumina-matrix.csv"               = "Illumina_matrix",
  "dna-methylation-scores.csv"        = "dnam_score",
  "camda-matrix.csv"                  = "camda_score",
  "estimate-cell-counts.csv"          = "ECC",
  "coverage-matrix.csv"               = "coverage_matrix",
  "cutadapt_filtered_reads_plot.txt"  = "reads"
)

files <- list()                              
files_path <- fread(file_list_path, header = FALSE)$V1

for (pth in files_path) {
  fname <- basename(pth)
  if (fname %in% names(key_map)) {
    key <- key_map[[fname]]
    files[[key]] <- pth
    assign(key, fread(pth) |> as.data.frame(), envir = .GlobalEnv)
    message(sprintf("Loaded %-30s → %s", fname, key))
  } else {
    message(sprintf("[skip] unmatched file: %s", fname))
  }
}

### Read & check phenotype ###
pheno <- data.frame(data.table::fread(phenotype.file))
expected_cols <- c("Sex","Country.x","Age_Rerc","Cancer_status","BMI_C",
                   "Smoke_status","Alc_Lifetime","Sample","batch")
missing <- setdiff(expected_cols, colnames(pheno))
if (length(missing) > 0) {
  stop("Phenotype file ", files$phenotype, 
       " is missing columns: ", paste(missing, collapse = ", "))
}
pheno$batch <- as.character(pheno$batch)

### Clean read counts ###
reads <- data.table::fread(files$reads)
setnames(reads, "Reads passing filters", "reads")
reads[, Sample := gsub("_R[12]$", "", Sample)]
reads <- unique(reads, by="Sample")
reads[, Lab_id := gsub(".*_(.+)", "\\1", Sample)]
reads[, Sample := gsub("-", ".", gsub("_.*","", Sample))]

pheno <- pheno %>% mutate(Sample=gsub("[_]+",".",Sample))
pheno <- merge(pheno, reads, by = "Sample", all.x = TRUE)
rownames(pheno) <- pheno$Sample

### Load cell counts and merge ###
ECC <- data.table::fread(files$ECC)
cell_types <- ECC$V1
ECC$V1 <- NULL
ECC <- t(ECC)
colnames(ECC) <- cell_types
rownames(ECC) <- gsub("_.*", "", rownames(ECC))
pheno <- cbind(pheno, ECC[pheno$Sample, , drop = FALSE])
pheno$Sample <- NULL

### Recode smoking ###
stopifnot("Smoke_status" %in% colnames(pheno))
stopifnot(all(pheno$Smoke_status %in% c("Former","Current","Never")))
pheno$Smoke_current_former <- as.integer(pheno$Smoke_status == "Current")
pheno$Smoke_current_former[pheno$Smoke_status=="Never"] <- NA
pheno$Smoke_former_never  <- as.integer(pheno$Smoke_status == "Former")
pheno$Smoke_former_never[pheno$Smoke_status=="Current"] <- NA
pheno$Smoke_current_never <- as.integer(pheno$Smoke_status == "Current")
pheno$Smoke_current_never[pheno$Smoke_status=="Former"] <- NA

### Read models ###
models <- read.csv(models.file, stringsAsFactors = FALSE)
models$model <- paste0(models$model, " + CD4T + CD8T + NK + Mono + Bcell + Granulocytes + reads")

### Helper: run ewaff ###
process_ewaff <- function(pheno, meth, out_folder, summary_folder=NULL, manifest=NULL, models) {
  colnames(meth) <- gsub("_.*", "", gsub("^X","",colnames(meth)))
  common <- intersect(rownames(pheno), colnames(meth))
  stopifnot(length(common) > 1)
  pheno <- pheno[common, , drop = FALSE]
  meth  <- as.matrix(meth[, common])
  results <- list()
  for (i in seq_len(nrow(models))) {
    var     <- models$var[i]
    mdl     <- models$model[i]
    nm      <- models$name[i]
    cat(Sys.time(), " — running model:", var, "\n")
    sites.ret <- ewaff::ewaff.sites(as.formula(mdl), variable.of.interest = var,
                                    methylation = meth, data = pheno, method = "glm")
    # summary + report if manifest provided
    if (!is.null(manifest)) {
      sum.ret <- ewaff::ewaff.summary(sites.ret, manifest$chr, manifest$start, meth)
      ewaff::ewaff.report(sum.ret, output.file = file.path(out_folder, paste0("report_", var, ".html")))
      results[[var]]$summary <- sum.ret
    }
    # write table
    if (!is.null(summary_folder)) dir.create(summary_folder, showWarnings=FALSE, recursive=TRUE)
    out_csv <- if (!is.null(summary_folder))
                   file.path(summary_folder, paste0("summary_statistics_", nm, ".csv"))
               else
                   file.path(out_folder, paste0("sites_ret_", nm, ".csv"))
    data.table::fwrite(cbind(features = rownames(sites.ret$table), sites.ret$table),
                       file = out_csv, row.names = FALSE)
    results[[nm]]$sites <- sites.ret
  }
  return(results)
}

### Helper: load matrices ###
prepare_data <- function(file, id_col=NULL, manifest_cols=NULL) {
  dt <- data.frame(data.table::fread(file))
  if (is.null(manifest_cols)) {
    manifest <- data.frame(chr="score", start=1:nrow(dt), end=1:nrow(dt))
  } else {
    manifest <- dt[, manifest_cols, drop = FALSE]
  }
  if (!is.null(id_col)) {
    rownames(dt) <- dt[[id_col]]
  } else {
    rownames(dt) <- paste(manifest$chr, manifest$start, sep=":")
  }
  list(data = dt, manifest = manifest)
}

### Run on each dataset ###
datasets <- list(
  Methylation_matrix = list(file = files$Methylation_matrix, manifest_cols = c("chr","start","end")),
  Illumina_matrix    = list(file = files$Illumina_matrix,    id_col="CpGs", manifest_cols = c("chr","start","end")),
  dnam_score         = list(file = files$dnam_score,         id_col="V1"),
  camda_score        = list(file = files$camda_score,        manifest_cols = c("chr","start","end"))
)

results <- lapply(names(datasets), function(name) {
  cat(Sys.time(), " — processing:", name, "\n")
  ds   <- datasets[[name]]
  outf <- file.path(output.dir, name)
  sumf <- file.path(outf, paste0("summary_statistics_", name))
  dir.create(outf, showWarnings=FALSE, recursive=TRUE)

  prep <- prepare_data(ds$file, ds$id_col, ds$manifest_cols)
  res  <- process_ewaff(pheno = pheno,
                        meth = prep$data,
                        out_folder = outf,
                        summary_folder = sumf,
                        manifest = prep$manifest,
                        models = models)

  save(res, file = file.path(outf, paste0("results_", name, ".rda")))
  return(res)
})

names(results) <- names(datasets)

