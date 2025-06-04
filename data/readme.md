## Target Regions: blood-cell-types-regions.bed

- This BED file contains the original **cell type-specific CpG marker regions**, directly downloaded from Loyfer et al., 2021 (Supplementary Table S4C). 
- No additional processing was performed on these regions.

- The `blood-cell-types-regions.bed` file is primarily used for **cell count estimation** (e.g., Houseman method) and related quality control (QC) steps in the pipeline.
- The file `blood_cell_types_extended.zip` is specifically generated for cell count estimation by processing Loyfer’s data and the regions in `blood-cell-types-regions.bed` using the [`make-blood-cell-bed.r`](https://github.com/MRCIEU/dnam-lung-cancer-pipeline/tree/main/scripts/make-blood-cell-bed.r) script.
- The script provides a reproducible way to generate the extended region set required for cell count estimation.

---

**Source:**  

Blood cell types extended list was downloaded from Loyfer et al., 2021 `Supplementary Table S4C`, which includes a list of 50,286 cell type-specific unmethylated markers.  
For this pipeline, the top 1,000 markers (hg19 coordinates) were used.
