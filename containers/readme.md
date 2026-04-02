# Container Definition Files

This folder contains Apptainer/Singularity definition files for DNA methylation and methyl-seq analysis environments.

---

## 1. meth-analysis/meth-analysis.def

- **Base image:** ubuntu:24.04
- **Purpose:** DNA methylation analysis with Python 3.12.3, Quarto 1.7.32 and R 4.5.3
- **Includes:** R packages (`data.table`, `ggplot2`, Bioconductor annotations), 
  Python `pandas` and `numpy`, GitHub R packages `meffonym` and `ewaff`.
  
The container is described as a sequence of incremental builds to be modularize the code
and reduce build time for development, especially debugging.

---

## 2. wgbs-analysis.def

- **Base image:** condaforge/mambaforge:latest  
- **Purpose:** Tools for methyl-seq analysis.  
- **Includes:** Genomics tools installed via mamba (e.g. Bwa-meth, Bowtie2, Picard, Methyldackel, Sambamba, Samtools).

---

## 3. qc-container.def

- **Base image:** rocker/r-ver:4.4.1  
- **Purpose:** R 4.4.1 environment with Quarto and essential data analysis packages.  
- **Includes:** Quarto (v1.6.42), core R packages (`data.table`, `ggplot2`, `dplyr`, etc.) for data processing and reporting.

---

## 4. wgbstools.def

- **Base image:** python:3.11-bullseye 
- **Purpose:** provides installation of wgbstools (https://github.com/nloyfer/wgbs_tools)
- **Includes:** wgbstools and Python 3.11 

(see [../data/readme.md](../data/readme.md) for how it can be used to extract 
cell-type specific DNA methylation described in Loyfer et al. 2023 https://doi.org/10.1038/s41586-022-05580-6)

---

## Notes

- Each container is designed for specific analysis steps and aims to ensure reproducibility and portability.  
- Build container images with:  
  
```bash
snakemake meth-analysis/meth-analysis.sif
snakemake wgbs-analysis.sif
snakemake qc-analysis.sif
snakemake wgbs-tools.sif
```

- Upload container images to docker hub:

```bash
# 1. create an account here: https://hub.docker.com/
# 2. login at the command-line
apptainer registry login --username <username> docker://docker.io
# 3.
apptainer push <container>.sif oras://docker.io/<username>/<container>:latest
```
