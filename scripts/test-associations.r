#!/usr/bin/env R

library(data.table)
library(ewaff)

args <- commandArgs(trailingOnly = TRUE)

file.list <- args[1]
phenotype.file  <- args[2]
models.file     <- args[3]
output.dir      <- args[4]

dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(output.dir, "sessionInfo.txt")); sessionInfo(); sink()

file.paths <- fread(file.list, header = FALSE)[[1]]

meth.file = file.paths[grep("methylation-matrix",file.paths)]
illumina.file = file.paths[grep("illumina-matrix",file.paths)]
dnam.scores.file = file.paths[grep("dna-methylation-scores",file.paths)]
camda.file = file.paths[grep("camda-matrix",file.paths)]
counts.file = file.paths[grep("estimate-cell-counts",file.paths)]
coverage.file = file.paths[grep("coverage-matrix",file.paths)]
reads.file = file.paths[grep("cutadapt_filtered_reads_plot",file.paths)]

## phenotypes ##
pheno <- data.frame(data.table::fread(phenotype.file))
stopifnot("Sample" %in% colnames(pheno))

## read counts ##
reads <- data.table::fread(reads.file)
setnames(reads, "Reads passing filters", "reads")
reads$Sample <- sub("_[^_]+$", "", reads$Sample)
reads = unique(reads)

pheno <- merge(pheno, reads, by = "Sample")
rownames(pheno) <- pheno$Sample

## cell counts ##
counts <- read.csv(counts.file,row.names=1)
counts <- t(counts)
pheno <- cbind(pheno, counts[pheno$Sample, , drop = FALSE])

pheno$Sample <- NULL

## prep datasets
prep <- function(x) {
  if ("chr" %in% colnames(x)) {
    manifest <- x[,c("chr","start","end")]
    dat <- as.matrix(x[,setdiff(colnames(x),colnames(manifest))])
    if (is.null(rownames(dat)))
      rownames(dat) <- with(manifest,paste0(chr,":",start,"-",end))
    list(manifest=manifest,dat=dat)
  } else {
    list(manifest=NULL,dat=as.matrix(x))
  }
}  

meth = prep(fread(meth.file,data.table=F))

illumina <- fread(illumina.file,data.table=F)
rownames(illumina) <- illumina$cpg
illumina$cpg <- NULL
illumina <- prep(illumina)

camda = prep(fread(camda.file,data.table=F))

dnam_scores = fread(dnam.scores.file,data.table=F)
rownames(dnam_scores) <- dnam_scores$name
dnam_scores$name <- NULL
dnam_scores <- prep(dnam_scores)

## models ##
models <- read.csv(models.file, stringsAsFactors = FALSE)

var.counts <- apply(counts,2,var,na.rm=T)
if (!all(is.na(var.counts))) {
  if (any(is.na(var.counts))) 
    cell.types <- colnames(counts)[!is.na(var.counts)]
  else
    cell.types <- colnames(counts)[-which.min(var.counts)]
} else
  cell.types <- c()

additional <- c("reads", cell.types)

models$model <- paste0(models$model, "+", paste(additional, collapse="+"))


## test associations
test_assocs <- function(pheno, dat, models, out.dir) {
  dir.create(out.dir, recursive=T)
  meth <- dat$dat
  manifest <- dat$manifest
  common <- intersect(rownames(pheno), colnames(meth))
  stopifnot(length(common) > 1)
  pheno <- pheno[match(common,rownames(pheno)), , drop=F]
  meth  <- as.matrix(meth[, match(common,colnames(meth)), drop=F])
  for (name in models$name) {
    var <- models$var[models$name == name]
    model <- models$model[models$name == name]
    cat(date(), name, var, "\n")
    stats = ewaff.sites(
      as.formula(model),
      variable.of.interest=var,
      methylation=meth,
      data=pheno,
      method="glm")
    stats$table$feature <- rownames(stats$table)
    if (!is.null(manifest)) {
      ret <- ewaff.summary(stats, manifest$chr,manifest$start,meth)
      ewaff.report(
        ret,
        output.file=file.path(out.dir,paste0("report-",name,".html")))
      stats$table = cbind(manifest,stats$table)                   
    }
    fwrite(stats$table,file=file.path(out.dir,paste0("stats-",name,".csv")), row.names=F)
    
    save(stats, file=file.path(out.dir, paste0("stats-",name,".rda")))
  }
}

test_assocs(pheno, meth, models, file.path(output.dir, "methylation"))
test_assocs(pheno, illumina, models, file.path(output.dir, "illumina"))
test_assocs(pheno, camda, models, file.path(output.dir, "camda"))
test_assocs(pheno, dnam_scores, models, file.path(output.dir, "dnam_scores"))

