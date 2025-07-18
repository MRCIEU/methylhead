#!/usr/bin/R

args <- commandArgs(trailingOnly = TRUE)
meth_file <- args[1]
output_file <- args[2]

assembly <- ifelse(length(args) > 2, args[3], "hg19")

stopifnot(assembly %in% c("hg19","hg38"))

if (assembly == "hg19") {
    annotations = c(
        "IlluminaHumanMethylation450kanno.ilmn12.hg19",
        "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
} else {
    annotations = c(
        "IlluminaHumanMethylationEPICv2anno.20a1.hg38")
}

for (annotation in annotations)
    library(annotation, character.only=T)
library(data.table)

meth_df <- data.frame(fread(meth_file))

annot <- lapply(annotations, function(annotation) getAnnotation(get(annotation)))
cols <- colnames(annot[[1]])
for (i in 1:length(annot)) cols <- intersect(cols, colnames(annot[[i]]))
annot <- lapply(annot, function(annot) annot[,cols])
annot <- do.call(rbind, annot)
annot$Name <- sub("_.*", "", annot$Name)
sites <- unique(annot$Name)
annot <- annot[match(sites, annot$Name),]

annot$loc <- paste(annot$chr, annot$pos)    
meth_df$loc <- paste(meth_df$chr, meth_df$end)
merged_data <- merge(meth_df, annot, by = c("loc"))
colnames(merged_data)[colnames(merged_data) == "chr.x"] <- "chr"
sample.ids<- colnames(meth_df)[!colnames(meth_df) %in% c("chr", "start", "end", "loc")]
selected_cols <- c("Name","chr","start","end",sample.ids)
methylation <- merged_data[, selected_cols]
Illumina_matrix <- methylation[, names(methylation) != "chr.y"]
colnames(Illumina_matrix)[1] <- "CpGs"
colnames(Illumina_matrix) <- gsub("\\.", "-", gsub("^X", "", colnames(Illumina_matrix)))
write.csv(Illumina_matrix, file = output_file,row.names=FALSE)
