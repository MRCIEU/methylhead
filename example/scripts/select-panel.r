#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
    stop("Usage: Rscript select-panel.r original-panel.csv read-counts.bed panel.csv panel.bed")
}

## original lung cancer target panel 
original_panel_filename = args[1]
## cell type regions with read counts
cell_types_filename = args[2]
## target panel csv for dataset 
panel_filename = args[3]
## target panel bed for dataset 
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

set.seed(1881)

selection = sample(1:nrow(targets), size=floor(nrow(targets)*0.8))

write.csv(
    targets[,selection], ## randomly remove 20%
    file=panel_filename,
    row.names=F, quote=F)

write.table(
    targets,
    file=bed_filename,
    sep="\t",
    row.names=F, col.names=F, quote=F)        

