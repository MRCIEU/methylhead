# Container Definition Files

This folder contains Apptainer/Singularity definition files for DNA methylation and WGBS analysis environments.

---

## 1. meth_analysis.def

- **Base image:** rocker/r-ver:4.4.3  
- **Purpose:** DNA methylation analysis with R 4.4.3 and Python integration.  
- **Includes:** R packages (`methylKit`, Bioconductor annotations), Python 3 with `pandas` and `numpy`, GitHub packages `meffonym` and `ewaff`.

---

## 2. wgbs_analysis.def

- **Base image:** condaforge/mambaforge:latest  
- **Purpose:** Tools for Whole-Genome Bisulfite Sequencing (WGBS) analysis.  
- **Includes:** Genomics tools installed via mamba (e.g. Bismark, Bowtie2, Picard, Methyldackel, Sambamba, Samtools).

---

## 3. quarto_analysis.def

- **Base image:** rocker/r-ver:4.4.1  
- **Purpose:** R 4.4.1 environment with Quarto and essential data analysis packages.  
- **Includes:** Quarto (v1.6.42), core R packages (`data.table`, `ggplot2`, `dplyr`, etc.) for data processing and reporting.

---

## Notes

- Each container is designed for specific analysis steps and aims to ensure reproducibility and portability.  
- Build containers with:  
  
```bash
  apptainer build <container>.sif <container>.def
```