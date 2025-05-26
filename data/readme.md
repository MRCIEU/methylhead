
## blood_cell_types_extended.bed

**How this file was built:**

For a full, reproducible recipe, see [`make-blood-cell-bed.r`](https://github.com/MRCIEU/dnam-lung-cancer-pipeline/tree/main/scripts/make-blood-cell-bed.r). This script walks through:

- Downloading and filtering the panel manifest (`panel-reduced.csv`)
- Extracting WGBS signal from public blood datasets using [`wgbs_tools`](https://github.com/nloyfer/wgbs_tools)
- Formatting to BED + β‑matrix for use with the Houseman cell count algorithm

You can use the same workflow as a template to generate reference files for other tissues or custom panels.

**Source:**  
Blood cell types extended list was downloaded from Loyfer et al., 2021 ([Supplementary Table S4C](https://www.nature.com/articles/s41586-021-03532-9)), which includes a list of 50,286 cell type-specific unmethylated markers.  
For this pipeline, the top 1,000 markers (hg19 coordinates) were used.
