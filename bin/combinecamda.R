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

processed_files <- c()  
sample.ids <- c()       
ci_list <- list()  

for (file in file_list) {
  sample_id <- gsub("_CpG_CAMDA.tsv", "", basename(file))
  data <- read.delim(file, header = TRUE)

  # CI verilerini olusturuyoruz ve ci_list'e ekliyoruz
  ci_data <- data %>%
    mutate(loc = paste(chr, pos, sep = "_")) %>%
    dplyr::select(loc, CI_lower, CI_upper) %>%
    mutate(sample_id = sample_id)
  ci_list[[sample_id]] <- ci_data  
  methylkit_data <- data %>%
    mutate(
      chrBase = paste0(chr, ".", pos),
      coverage = rep(11, nrow(data)), 
      freqC = ratio * 100,  
      freqT = 100 - freqC,  
      strand = ifelse(strand == "+", "F", "R"),
      loc = paste(chr, pos, sep = "_")
    ) %>%
    dplyr::select(chrBase, chr, pos, strand, coverage, freqC, freqT, loc)  

  # Islenmis dosyayi geçici bir dizine kaydediyoruz
  processed_file <- paste0(tempdir(), "/", sample_id, "_processed.tsv")
  write.table(methylkit_data, processed_file, sep = "\t", row.names = FALSE, quote = FALSE) 
  processed_files <- c(processed_files, processed_file)
  sample.ids <- c(sample.ids, sample_id)
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
