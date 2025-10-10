#!/usr/bin/env Rscript

args         <- commandArgs(trailingOnly = TRUE)
meth_file    <- args[1]
output_file  <- args[2]
sites_file   <- args[3]

library(data.table)
library(meffonym)

meth <- fread(meth_file,data.table=F)
rownames(meth) <- meth$cpg
meth <- meth[,setdiff(colnames(meth), c("cpg","chr","start","end"))]
meth <- as.matrix(meth)

models <- meffonym.models()
outputs <- sapply(meffonym.models(), function(model) {
  try(meffonym.score(meth,model), silent=T)
}, simplify=F)
is.good <- sapply(outputs, length) > 1
scores <- do.call(rbind, lapply(outputs[is.good],function(out) out$score))
colnames(scores) <- colnames(meth)
sites <- data.frame(
  name=rownames(scores),
  sites=sapply(outputs[is.good],function(out) length(out$vars)),
  used=sapply(outputs[is.good],function(out) length(out$sites)))

scores <- data.frame(name=rownames(scores), scores, check.names=F)
write.csv(scores,file=output_file,row.names=F)
write.csv(sites,file=sites_file,row.names=F)
