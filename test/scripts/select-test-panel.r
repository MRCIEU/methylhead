#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
    stop("Usage: Rscript select-test-panel.r original-panel.csv read-counts.bed panel.csv panel.bed")
}

## original lung cancer target panel (hg19 coords)
original_panel_filename = args[1]
## cell type regions with read counts (hg19 coords)
cell_types_filename = args[2]
## target panel csv for test dataset (new coords and hg19 coords)
panel_filename = args[3]
## target panel bed for test dataset (hg19 coords)
bed_filename = args[4]            

options(scipen=999)

## load cell count regions with counts info
## (keep only regions with at least 5bp coverage)
cell_type_regions = read.table(cell_types_filename, sep="\t", header=F)
read_counts = cell_type_regions[,-(1:3)]
cell_type_regions = cell_type_regions[,1:3]
colnames(cell_type_regions) = c("chr","start","end")
colnames(read_counts) = paste0("s",1:ncol(read_counts))
coverage = read_counts/(cell_type_regions$end-cell_type_regions$start)
mincov = apply(coverage, 1, min)
cell_type_regions = cell_type_regions[mincov > 5,]

## load original lung cancer panel
original_panel = read.csv(original_panel_filename)

## merge panel and filtered cell count regions
targets = rbind(cell_type_regions, original_panel[,c("chr","start","end")])
targets = unique(targets)

## extend target regions by 500 bp
targets$start = targets$start - 500
targets$end = targets$end + 500

stopifnot(all(targets$start > 0))

## merge overlapping regions
targets = targets[order(targets$chr, targets$start),]
targets$keep = F
targets$old_start = targets$start
targets$old_end = targets$end
chr = ""
start = end = NA
for (i in 1:nrow(targets)) {
    if (targets$chr[i] == chr && targets$start[i] <= end)
        end = max(end, targets$end[i])
    else {
        if (i > 1) {
            targets$keep[i-1] = TRUE
            targets$start[i-1] = start
            targets$end[i-1] = end
        }
        chr = targets$chr[i]
        start = targets$start[i]
        end = targets$end[i]
    }
}
targets$keep[nrow(targets)] = TRUE
targets$start[nrow(targets)] = start
targets$end[nrow(targets)] = end
targets = targets[targets$keep,]
targets = targets[,c("chr","start","end")]

## convert targets to new genome coordinates
converted_targets = targets
converted_targets$original_start = targets$start
converted_targets$original_end = targets$end
for (chr in unique(targets$chr)) {
    idx = which(targets$chr == chr)
    target_length = targets$end[idx]-targets$start[idx]
    converted_targets$end[idx] = cumsum(target_length)
    converted_targets$start[idx] = converted_targets$end[idx]-target_length
}

set.seed(1881)

# The following block of code creates a smaller test panel by taking a random 80% sample from the full set of target regions. 
# This is done to ensure the final test dataset is a manageable size while still being representative of the original data.

# Get the total number of rows from the full targets data frame.
total_rows <- nrow(converted_targets)

#Calculate the desired sample size (80%) and round down.
sample_size <- floor(total_rows * 0.80)

# Generate a random sample of row indices. Because we used set.seed(), this will produce the same indices every time.
random_indices <- sample(1:total_rows, size = sample_size)

# Create the new, smaller panel_targets data frame using the same reproducible indices.
panel_targets <- converted_targets[random_indices, ]

## save panel
write.csv(
    panel_targets,
    file=panel_filename,
    row.names=F, quote=F)

write.table(
    converted_targets[,c("chr","original_start","original_end")],
    file=bed_filename,
    sep="\t",
    row.names=F, col.names=F, quote=F)        

