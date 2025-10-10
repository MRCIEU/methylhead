# methylhead · Panel‑WGMS Analysis Pipeline

<p align="center">
  <img src="flowchart/methylhead.-logo.png" alt="methylhead logo" width="220">
</p>

[![CI](https://github.com/MRCIEU/methylhead/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/MRCIEU/methylhead/actions/workflows/ci.yml)
[![Nextflow ≥23.10.0](https://img.shields.io/badge/nextflow-%E2%89%A523.10.0-brightgreen)](https://www.nextflow.io/)
[![Docker Pulls](https://img.shields.io/docker/pulls/onuroztornaci/methylhead-pipeline)](https://hub.docker.com/r/onuroztornaci/methylhead-pipeline)
![Apptainer](https://img.shields.io/badge/apptainer-SIF-blue)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**methylhead** is a modular **nextflow** workflow that turns raw targeted methyl-seq FASTQ files into QC‑checked methylation matrices, cell‑composition estimates and model‑based dnam risk scores—ready for statistics or reporting.

---
## 🌟 Why methylhead? — Feature highlights

| Feature                   |  Description                                                       |
|---------------------------|--------------------------------------------------------------------|
| End‑to‑end panel‑WGBS     |  From raw FASTQ to sample‑level dnam risk scores with a single command |
| Cell‑composition inference|  Blood‑cell deconvolution using bundled reference libraries       |
| Model‑based predictions   |  Runs arbitrary EWAS/age/risk models defined in a CSV            |
| Reproducible & portable   |  Fully containerised (Apptainer); no system installation |
| Modular Nextflow core     |  Parallel execution, --resume, profile support                   |
| Rich QC out‑of‑the‑box    |  Per‑sample & per‑locus thresholds, MultiQC and Quarto HTML/PDF reports |

---
---

## · Prerequisites

| Requirement   | Tested version | Check with            |
| ------------- | -------------- | --------------------- |
| **Apptainer** |  ≥ 1.1.0       | `apptainer --version` |
| **Conda**     |  ≥ 23.x        | `conda -V`            |
| Internet      | outbound HTTPS | —                     |

* **Apptainer ≥ 1.1** ([install guide](https://apptainer.org/docs/))
* **Conda ≥ 23.x** ([install guide](https://docs.conda.io/en/latest/miniconda.html))

---

## 1 · Clone the repository

```bash
git clone git@github.com:MRCIEU/methylhead.git
cd methylhead
```

---

## 2 · (One‑off) Build the reference genome (≈ 2 h)

```bash
bash scripts/create-reference.sh hg19 PATH/TO/HG19/GENOME/INDEX
```

* `PATH/TO/HG19/GENOME/INDEX` The location where the human genome (hg19) should be downloaded and indexed.

You can skip this step if you already have an indexed hg19 reference.

---

## 3 · Install nextflow using conda

If nextflow is not already installed, you can install it using conda.

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict
conda create -y -n methylhead nextflow -c bioconda
conda activate methylhead
```

## 4 · Quick start with test dataset

A [small test dataset](test/readme.md) has been created in the [test/](test) directory
to quickly verify that your environment is ready.


## 5 · Run on your own samples

Please first edit as appropriate directory paths and settings in the first section of 
[nextflow.config](nextflow.config). 
See this [link](https://www.nextflow.io/docs/latest/config.html) for further information 
about the configuration file.

```bash
nextflow -C nextflow-test.config run main.nf -resume
```

* [-resume](https://nextflow.io/docs/latest/cache-and-resume.html) lets nextflow 
  **pick up from where a previous run left off**

## 5 ·  Mandatory inputs

| ---------------- | ---------------------------------------------------------- | ----------------------------- |
| `samplesheet`   | list of paired fastq files for each sample                   | `samplesheet.csv`    |
| `genome_fasta`  | genome sequencing indexed for alignment                      | `genome/hg19.fa`                    |
| `cell_reference`| cell-type reference dataset for estimating cell counts        | `data/blood-cell-type-reference.csv.gz`|
| `panel`         | genomic target regions                                       | `panel.csv`                            |
| `phenotype`     | sample‑level metadata                                   | `phenotype.csv`                       |
| `models`        | models for association testing                               | `models.csv`                      |
| `assembly`      | genomic assembly (must be hg19 or hg38)                      | `hg19`                    |

> **See [`input/readme.md`](input/readme.md) for file formats & examples.**   

---

## 6 · Workflow overview

This directory contains a single file:

| File           | Description                   |
| -------------- | ---------------------------- |
| workflow.png   | Auto-generated Nextflow DAG   |

The [`workflow.png`](flowchart/workflow.png) file visualizes the task-level dependencies in the pipeline, as produced by `nextflow dag`.
> **See [`/flowchart/readme.md`](flowchart/readme.md) for file formats step by step.**
---

## 7 · Containers in use

| Flag         | Default URI                                                        | Includes                        |
| ------------ | ------------------------------------------------------------------ | ------------------------------- |
| `wgbs_image` | `oras://docker.io/onuroztornaci/methylhead-pipeline:wgbs_analysis` | WGBS aligners & QC              |
| `meth_image` | `oras://docker.io/onuroztornaci/methylhead-pipeline:meth_analysis` | R 4.4.3, Python 3, Bioconductor |
| `qc_image`   | `oras://docker.io/onuroztornaci/methylhead-pipeline:qc_container`  | R 4.4.1, Quarto                 |

Build your own images → see [`/container-def-files`](/container-def-files/readme.md).

---

## 8 · Utilities

* [methylkit-to-matrix](scripts/methylkit-to-matrix/) Script for extracting a basic dataset from [MethylKit](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-10-r87) output files

---

## 9 · Troubleshooting cheatsheet

| Symptom                       | Likely cause & fix                                                        |
| ----------------------------- | ------------------------------------------------------------------------- |
| `ERROR: Apptainer not found`  | Install Apptainer ≥ 1.1 and add it to `$PATH`.                            |
| Java `<11` warning            | Forgot to `conda activate methylhead`.                                    |
| `Index not found for hg19.fa` | Run **4 · reference build** or point `--genome_fasta` to an indexed ref.  |
---

## 10 · Pipeline Status
 
The badges above provide a real-time summary of the pipeline's continuous integration status, software requirements, container support, and license.

Happy methylating 🧬🚀
