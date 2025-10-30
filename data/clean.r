#!/usr/bin/env Rscript

library(data.table)

args <- commandArgs(trailingOnly = TRUE)

ref_file <- args[1] ## blood cell type reference data file
output_file <- args[2] ## cleaned data

ref_df <- fread(ref_file, data.table=F)

regions <- ref_df[,c("chr","start","end"), drop=F]

ref_df <- ref_df[, setdiff(names(ref_df), c("chr","start","end","startCpG","endCpG")), drop=F]

clean_names <- names(ref_df)

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
clean_names <- sapply(clean_parts, paste, collapse=".")

# Collapse replicates by group (rowMeans if multiple replicates)
cols_by_group <- split(seq_along(clean_names), clean_names)
clean_df <- sapply(cols_by_group, function(idx) {
  if (length(idx) == 1) ref_df[[idx]]
  else rowMeans(as.matrix(ref_df[, idx, drop=FALSE]))
})
clean_df <- cbind(regions,clean_df)

fwrite(clean_df, file=output_file, row.names=F)
