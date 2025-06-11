#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

input_panel_filename = args[1]
coverage_filename = args[2]
panel_filename = args[3]
bed_filename = args[4]

## load cell count regions with coverage info
## (keep only regions with at least 5bp coverage)
coverage = read.table(coverage_filename, sep="\t", header=F)
cell_type_regions = coverage[,1:3]
coverage = coverage[,-(1:3)]
colnames(cell_type_regions) = c("chr","start","end")
colnames(coverage) = paste0("s",1:ncol(coverage))
coverage = coverage/(cell_type_regions$end-cell_type_regions$start)
mincov = apply(coverage, 1, min)
cell_type_regions = cell_type_regions[mincov > 5,]

## load original lung cancer panel
panel = read.csv(input_panel_filename)
panel$chr = paste0("chr",panel$chr)

## merge panel and filtered cell count regions
targets = rbind(cell_type_regions, panel[,c("chr","start","end")])
targets = unique(targets)

## extend target regions by 100 bp
targets$start = targets$start - 100
targets$end = targets$end + 100

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
    target_length = targets$end[idx]-targets$start[idx] + 1
    converted_targets$end[idx] = cumsum(target_length)
    converted_targets$start[idx] = converted_targets$end[idx] - target_length + 1
}

## save panel
write.csv(
    converted_targets,
    file=panel_filename,
    row.names=F, quote=F)

write.table(
    converted_targets[,c("chr","original_start","original_end")],
    file=bed_filename,
    sep="\t",
    row.names=F, col.names=F, quote=F)        

