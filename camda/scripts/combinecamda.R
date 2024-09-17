library(dplyr)
library(methylKit)

args <- commandArgs(trailingOnly = TRUE)
samples <- args[1]
output_file <- args[2]

file <-read.csv(samples,header=F)
file<-data.frame(file)
file_paths <- file$V2
cleaned_paths <- gsub("\\[|\\]", "", file_paths)
file_list <- trimws(cleaned_paths)

processed_files <- c()  
sample.ids <- c()       
for (file in file_list) {
sample_id <- gsub("_CpG_CAMDA.tsv", "", basename(file))
 data <- read.delim(file, header = TRUE) 
      data <- data %>%
        filter((CI_upper - CI_lower) <  0.40) 
  methylkit_data <- data %>%
    mutate(
      chrBase = paste0(chr, ".", pos),
      coverage = rep(11, nrow(data)), 
      freqC = ratio * 100,  
      freqT = 100 - freqC,  
      strand = ifelse(strand == "+", "F", "R")
    ) %>%
    dplyr::select(chrBase, chr, pos, strand, coverage, freqC, freqT)
  processed_file <- paste0(tempdir(), "/", sample_id, "_processed.tsv")
  write.table(methylkit_data, processed_file, sep = "\t", row.names = FALSE, quote = FALSE)
  processed_files <- c(processed_files, processed_file)
  sample.ids <- c(sample.ids, sample_id)
}
file.list<-as.list(processed_files)
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
 dplyr::select(-matches("coverage|numCs|numTs")) 
write.csv(camda_matrix, output_file, row.names = FALSE)
