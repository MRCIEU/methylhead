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
meth.seqlm.file   <- args[7]
counts.file     <- args[8]
coverage.file   <- args[9]
reads.file      <- args[10]
output.dir      <- args[11]

dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(output.dir, "sessionInfo.txt")); sessionInfo(); sink()

## phenotypes ##
pheno <- data.frame(data.table::fread(phenotype.file))
stopifnot("sample_id" %in% colnames(pheno))

## read counts ##
reads <- data.table::fread(reads.file)
setnames(reads, "Reads passing filters", "reads")
reads$sample_id <- sub("_R[12]+$", "", reads$Sample)
reads$Sample <- NULL
reads = unique(reads)

pheno <- merge(pheno, reads, by = "sample_id")

## cell counts ##
counts <- read.csv(counts.file,row.names=1,check.names=F)
counts <- t(counts)
pheno <- cbind(
  pheno,
  counts[match(pheno$sample_id,rownames(counts)),])

rownames(pheno) <- pheno$sample_id
pheno$sample_id <- NULL

## prep datasets
prep <- function(x,ids,row.names=NULL) {
  if (!is.null(row.names))
    rownames(x) <- x[[row.names]]
  if ("chr" %in% colnames(x)) {
    manifest <- x[,c("chr","start","end")]
    if (is.null(row.names)) 
      rownames(x) <- with(manifest,paste0(chr,":",start,"-",end))
  } else 
    manifest <- data.frame(chr=rownames(x),start=1,end=1)
  dat <- as.matrix(x[,intersect(colnames(x),ids)])
  rownames(dat) <- rownames(manifest) <- rownames(x)
  list(manifest=manifest,dat=dat)
}  

meth <- prep(fread(meth.file,data.table=F),rownames(pheno))

illumina <- prep(fread(illumina.file,data.table=F),rownames(pheno),"cpg")

camda <- prep(fread(camda.file,data.table=F),rownames(pheno))

meth.seqlm <- prep(fread(meth.seqlm.file,data.table=F),rownames(pheno))

dnam.scores <- prep(fread(dnam.scores.file,data.table=F),rownames(pheno),"name")

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
    if (!grepl("^methylation", model)
        && length(unique(na.omit(pheno[[var]]))) == 2) {
      family <- "binomial"
      if (is.character(pheno[[var]]))
        pheno[[var]] <- as.factor(pheno[[var]])
    } else
      family <- "gaussian"
    random.subset = ifelse(0.05*nrow(meth) > 3, 0.05, 1)
    stats = ewaff.sites(
      as.formula(model),
      variable.of.interest=var,
      methylation=meth,
      data=pheno,
      method="glm",
      family=family,
      random.subset=random.subset)
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
test_assocs(pheno, dnam.scores, models, file.path(output.dir, "dnam.scores"))
test_assocs(pheno, meth.seqlm, models, file.path(output.dir, "methylseqlm"))


