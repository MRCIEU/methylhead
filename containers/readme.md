# Container Definition Files

This folder contains Apptainer/Singularity definition files for DNA methylation and methyl-seq analysis environments.

---

## 1. meth-analysis.def

- **Base image:** rocker/r-ver:4.4.3  
- **Purpose:** DNA methylation analysis with R 4.4.3 and Python integration.  
- **Includes:** R packages (`methylKit`, Bioconductor annotations), Python 3 with `pandas` and `numpy`, GitHub packages `meffonym` and `ewaff`.

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
apptainer build <container>.sif <container>.def
```

- Upload container images to docker hub:

```bash
# 1. create an account here: https://hub.docker.com/
# 2. login at the command-line
apptainer registry login --username <username> docker://docker.io
# 3.
apptainer push <container>.sif oras://docker.io/<username>/<container>:latest
```
