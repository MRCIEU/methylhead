
### How `blood_cell_types_extended.bed` was built

For a full, reproducible recipe see **[`make-blood-cell-bed.r`](https://github.com/MRCIEU/dnam-lung-cancer-pipeline/tree/main/scripts/)**. It walks through:

1. Downloading & filtering the panel manifest (`panel-reduced.csv`).
2. Extracting WGBS signal from public blood datasets with [`wgbs_tools`](https://github.com/nloyfer/wgbs_tools).
3. Formatting to BED + β‑matrix for use with the Houseman cell‑count algorithm.

* You can use the same workflow as a template to generate reference files for other tissues or custom panels.

- Blood cell types extended download from [loyfer's paper](https://www.nature.com/articles/s41586-022-05580-6): Supplementary Tables 1–17. Supplementary Table S4C.
- List of 50286 cell type-specific unmethylated markers (top 1000, hg19) was used.





