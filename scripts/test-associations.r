#!/usr/bin/env R

library(data.table)
library(ewaff)

args <- commandArgs(trailingOnly = TRUE)

phenotype.file  <- args[1]
models.file     <- args[2]
meth.file       <- args[3]
illumina.file   <- args[4]
dnam.scores.file <- args[5]
camda.file      <- args[6]
counts.file     <- args[7]
coverage.file   <- args[8]
reads.file      <- args[9]
output.dir      <- args[10]

dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(output.dir, "sessionInfo.txt")); sessionInfo(); sink()

## phenotypes ##
pheno <- data.frame(data.table::fread(phenotype.file))
stopifnot("sample_id" %in% colnames(pheno))

## read counts ##
reads <- data.table::fread(reads.file)
setnames(reads, "Reads passing filters", "reads")
reads$sample_id = reads$Sample

idx = lapply(pheno$sample_id, grep, x=reads$sample_id)
idx = cbind(
  pheno=rep(1:nrow(pheno),sapply(idx,length)),
  reads=unlist(idx))
reads = data.frame(
  sample_id=pheno$sample_id[idx[,"pheno"]],
  reads=reads$reads[idx[,"reads"]])

reads = unique(reads)

pheno <- merge(pheno, reads, by = "sample_id")
rownames(pheno) <- pheno$sample_id

## cell counts ##
counts <- read.csv(counts.file,row.names=1,check.names=F)
counts <- t(counts)
pheno <- cbind(pheno, counts[pheno$sample_id, , drop = FALSE])

pheno$sample_id <- NULL

## prep datasets
prep <- function(x) {
  if ("chr" %in% colnames(x)) {
    manifest <- x[,c("chr","start","end")]
    dat <- as.matrix(x[,setdiff(colnames(x),colnames(manifest))])
    if (is.null(rownames(dat)))
      rownames(dat) <- with(manifest,paste0(chr,":",start,"-",end))
    list(manifest=manifest,dat=dat)
  } else {
    manifest = data.frame(chr=rownames(x),start=1,end=1)
    list(manifest=manifest,dat=as.matrix(x))
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

