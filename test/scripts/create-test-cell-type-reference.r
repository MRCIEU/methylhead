#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#args = c("../data/blood-cell-type-reference.csv.gz", "input/panel.csv", "data/blood-cell-type-reference.csv")

original_reference_filename = args[1]
panel_filename = args[2]
cell_type_reference_filename = args[3]

panel = read.csv(panel_filename)

original_reference = read.table(original_reference_filename,header=T)
panel_loc = sapply(1:nrow(original_reference), function(i) {
    panel_loc = which(
        original_reference$chr[i] == panel$chr
        & original_reference$start[i] >= panel$original_start
        & original_reference$end[i] <= panel$original_end)
    if (length(panel_loc) == 0)
        NA
    else
        panel_loc[1]
})
converted_start = panel$start[panel_loc] + original_reference$start - panel$original_start[panel_loc]
converted_end = converted_start + original_reference$end - original_reference$start
converted_reference = original_reference
converted_reference$start = converted_start
converted_reference$end = converted_end

converted_reference = converted_reference[!is.na(panel_loc),]

write.csv(converted_reference, file=cell_type_reference_filename, row.names=F)
