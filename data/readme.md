# Blood cell-type reference dataset

This document outlines the process to reproduce the blood cell type reference
(`blood-cell-type-reference.csv.gz`). 
This file contains blood cell type-specific data 
in cell-type specific genomic regions (hg19 coordinates) listed in 
Supplementary Table S4C of [Loyfer et al., 2023](https://pubmed.ncbi.nlm.nih.gov/36599988/).

---

### Software 

## Steps

### 0. Install wgbstools

Either install `wgbstools` as described (https://github.com/nloyfer/wgbs_tools)
or pull an apptainer image with `wgbstools` already installed:

```
apptainer pull wgbstools.sif oras://docker.io/matthewsuderman/wgbstools
```

The following environment variable provides a convenient way to 
invoke wgbstools via apptainer:
```
GENOME_DIR=path/to/directory/of/hg19.fa
WGBSTOOLS="apptainer run -B $GENOME_DIR:/opt/wgbstools/references wgbstools.sif wgbstools"
```

### 1. Download Data

Download cell type-specific data from the Loyfer et al. dataset for all blood cell types.

```
bash ../scripts/blood-cell-type-download.sh blood-data
```

### 2. Select target regions

Extract the 'chr','start' and 'end' CpG sites from Table S4C from the Loyfer et al. dataset
(https://pmc.ncbi.nlm.nih.gov/articles/instance/9834054/bin/41586_2022_5580_MOESM4_ESM.xlsx). 
Save these as [cell-type-regions.bed](cell-type-regions.bed) 
after sorting by chr/start.

### 3. Index target regions 

Additional columns must be added to the bed file by `wgbstools` prior to data extraction.

```
$WGBSTOOLS convert -L cell-type-regions.bed --out_path cell-type-regions-wgbs.bed
```

### 4. Extract blood cell type-specific DNA methylation data

```
$WGBSTOOLS beta_to_table cell-type-regions-wgbs.bed --betas blood-data/*  \
	| column -t | gzip -c > blood-cell-type-reference-raw.csv.gz
```

### 5. Clean reference data

```
Rscript clean.r blood-cell-type-reference-raw.csv.gz blood-cell-type-reference.csv.gz
```

---

### Reference

* **Loyfer, N., … Kaplan, T. (2023).** A DNA methylation atlas of normal human cell types. *Nature, 613*(7943), 355–364. (PMID [36599988](https://pubmed.ncbi.nlm.nih.gov/36599988/)).
