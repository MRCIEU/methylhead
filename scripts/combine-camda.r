#!/usr/bin/env Rscript

library(dplyr)
library(methylKit)

args <- commandArgs(trailingOnly = TRUE)
samples <- args[1] 
output_file <- args[2]

file <- read.csv(samples, header = F)
file <- data.frame(file)
file_paths <- file$V2
cleaned_paths <- gsub("\\[|\\]", "", file_paths)
file_list <- trimws(cleaned_paths)

processed_files <- vector("character", length(file_list))
sample.ids <- vector("character", length(file_list))
ci_list <- vector("list", length(file_list))
names(ci_list) <- basename(file_list)

for (i in seq_along(file_list)) {
  file_path <- file_list[i]
  sample_id <- gsub("_CpG_CAMDA.tsv", "", basename(file_path))
  data <- data.table::fread(file_path)
  data[, loc := paste(chr, pos, sep = "_")]
  ci_data <- data[, .(loc, CI_lower, CI_upper)]
  ci_data[, sample_id := sample_id]
  ci_list[[sample_id]] <- ci_data
  methylkit_data <- data[, .(
    chrBase = paste0(chr, ".", pos),
    chr,
    pos,
    strand = ifelse(strand == "+", "F", "R"),
    coverage = 11L,
    freqC = ratio * 100,
    freqT = 100 - (ratio * 100),
    loc
  )]
  processed_file <- file.path(tempdir(), paste0(sample_id, "_processed.tsv"))
  data.table::fwrite(methylkit_data, processed_file, sep = "\t", quote = FALSE)
  processed_files[i] <- processed_file
  sample.ids[i] <- sample_id
}

file.list <- as.list(processed_files)
sample.ids <- as.list(sample.ids)
myobj <- methRead(file.list,
                  sample.id = sample.ids,
                  pipeline = "amp",                
                  assembly = "hg19",
                  treatment = rep(0, length(sample.ids)),
                  mincov = 10)

camda_matrix <- unite(myobj, destrand = FALSE, min.per.group = as.integer(length(sample.ids) * 0.50))

pm <- percMethylation(camda_matrix)
pm <- pm / 100
camda_matrix <- data.frame(camda_matrix)
camda_matrix <- cbind(camda_matrix, pm)

camda_matrix <- camda_matrix %>%
  mutate(loc = paste(chr, start, sep = "_")) 

camda_matrix <- camda_matrix %>%
  dplyr::select(-matches("coverage|numCs|numTs"))

for (sample_id in sample.ids) {
  ci_data <- ci_list[[sample_id]]
  camda_matrix <- camda_matrix %>%
    left_join(ci_data, by = "loc") %>%
    rename_with(~paste0(sample_id, "_", .), c("CI_lower", "CI_upper"))
}

camda_matrix <- camda_matrix %>%
  dplyr::select(-matches("sample_id|strand|loc"))
write.csv(camda_matrix, output_file, row.names = FALSE)
file.remove(processed_files)