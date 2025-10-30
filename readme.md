# Methylhead Pipeline - Snakemake Implementation

Snakemake conversion of the methylhead Nextflow pipeline for targeted methyl-seq analysis.

## Overview

This pipeline processes targeted methyl-seq data from raw FASTQs through to:
- Quality control and trimming
- DNA methylation-aware alignment (BWA-meth and BSMAP)
- Methylation calling
- Cell type deconvolution
- DNA methylation scores
- Association testing

## Requirements

- **Snakemake** ≥ 7.0
- **Python** ≥ 3.8 with packages **pandas** and **yaml**
- **Apptainer** ≥ 1.1.0

## Installation

```bash
# Install Snakemake via conda
conda create -n snakemake -c conda-forge -c bioconda snakemake pandas
conda activate snakemake

# Clone repository
git clone <repository-url>
```

You can verify correct installation by 
analysing the [test dataset](test). 

## Quick Start

### 1. Prepare inputs

Pipeline inputs include 
`samplesheet.csv`, `phenotypes.csv`, `models.csv`, and `panel.csv`.
See our [guide](input/readme.md) for preparing these inputs and examples. 

### 2. Configure the pipeline

Create a pipeline configuration file `config.yaml` 
by modifying [config-example.yml](config-example.yml).

The main modifications are in the 'inputs' and 'paths' sections.

### 3. Run the Pipeline

#### Stand-alone server

```bash
```

#### Cluster (SLURM)

```bash
```

## License

Same as original methylhead pipeline.
