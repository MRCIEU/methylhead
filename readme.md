# methylhead · Panel‑WGBS Analysis Pipeline
[![CI](https://github.com/MRCIEU/methylhead/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/MRCIEU/methylhead/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow-23.10.0-brightgreen)](https://www.nextflow.io/)
[![Docker Pulls](https://img.shields.io/docker/pulls/onuroztornaci/methylhead-pipeline)](https://hub.docker.com/r/onuroztornaci/methylhead-pipeline)
![Apptainer](https://img.shields.io/badge/apptainer-SIF-blue)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**methylhead** is a modular **Nextflow** workflow that turns raw targeted‑bisulfite FASTQ files into QC‑checked methylation matrices, cell‑composition estimates and model‑based risk scores—ready for statistics or reporting.

---
## 🌟 Why methylhead? — Feature highlights

| Feature                   |  Description                                                       |
|---------------------------|--------------------------------------------------------------------|
| End‑to‑end panel‑WGBS     |  From raw FASTQ to sample‑level risk scores with a single command |
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
# Pick any folder you like
git clone git@github.com:MRCIEU/methylhead.git
cd methylhead
```

---

## 2 · Quick start (≈ 5-10 min)

```bash
# Install & activate Nextflow if you haven’t yet
conda create -y -n methylhead nextflow -c bioconda
conda activate methylhead

# Run the built‑in demo (downloads containers on first run)
nextflow -C nextflow-test.config run main.nf 
```
* -C <file> tells Nextflow to **merge** the specified config file with the default nextflow.config. More: [Nextflow docs › configuration](https://www.nextflow.io/docs/latest/config.html)
* The demo dataset is documented inside the test/ folder—see [test/readme.md](test/readme.md) for details.


## 3 · (One‑off) Build the reference genome (≈ 2 h)

```bash
bash scripts/create-reference.sh -N you@example.com
```

*Creates `reference/hg19/` with all **bwameth** indices.*
Skip this step if you already have an indexed hg19 reference.

---

## 4 · Run on your own samples

```bash
nextflow run main.nf \
  --data            path/to/fastqs/*.fastq.gz \
  --genome_folder   path/to/hg19.fa \
  --cell_reference  path/to/cell-reference.csv \
  --panel           path/to/panel.csv \
  --phenotype       path/to/phenotype.csv \
  --models          path/to/models.csv \
  --outdir          results/ \
  -N you@example.com \
  --resume
```
* Leave out -N if you do **not** want an email summary.
* --resume lets Nextflow **pick up from where a previous run left off**—it will skip any steps that already finished successfully. More: [Nextflow docs › resume](https://nextflow.io/docs/latest/cache-and-resume.html)

### Mandatory parameters

| Flag              | Description                                                | Example                                |
| ----------------- | ---------------------------------------------------------- | -------------------------------------- |
| `--data`          | Glob of **gz‑compressed FASTQ** files                      | `mydata/*.fastq.gz`                    |
| `--genome_folder` | *Indexed* hg19 FASTA (`.fa` + `.bwt/.amb/...`)             | `reference/hg19.fa`                    |
| `--cell_reference`| cell-type-specific reference for cell-count estimation     | `data/blood-cell-type-reference.csv.gz`|
| `--panel   `      | CSV with per‑locus QC thresholds                           | `panel.csv`                            |
| `--phenotype`     | Sample‑level metadata                                      | `pheno.csv`                            |
| `--models`        | EWAS / risk‑prediction model definitions                   | `models.csv`                           |
> **See [`input/readme.md`](input/readme.md) for file formats & examples.**   

Optional flags:

| Flag                | Purpose                 | Default    |
| ------------------- | ----------------------- | ---------- |
| `--outdir`          | Where results go        | `results/` |
| `-N`                | Email run summary       | off        |
| `--wgbs_image` etc. | Override container URIs | built‑ins  |

---

## 5 · Outputs at a glance

```
results/
├── alignments/          # deduplicated BAM + stats
├── methylation_calls/   # BedGraphs per sample
├── matrices/            # CpG, coverage & 450k matrices
├── qc/                  # MultiQC + HTML/PDF report
└── predictions/         # Risk scores & association tests
```

---

## 6 · Workflow overview

This directory contains a single file:

| File           | Description                   |
| -------------- | ---------------------------- |
| workflow.png   | Auto-generated Nextflow DAG   |

The [`workflow.png`](/flowchart/workflow.png) file visualizes the task-level dependencies in the pipeline, as produced by `nextflow dag`.
> **See [`/flowchart/readme.md`](/flowchart/readme.md) for file formats step by step.**
---

## 7 · Containers in use

| Flag         | Default URI                                                        | Includes                        |
| ------------ | ------------------------------------------------------------------ | ------------------------------- |
| `wgbs_image` | `oras://docker.io/onuroztornaci/methylhead-pipeline:wgbs_analysis` | WGBS aligners & QC              |
| `meth_image` | `oras://docker.io/onuroztornaci/methylhead-pipeline:meth_analysis` | R 4.4.3, Python 3, Bioconductor |
| `qc_image`   | `oras://docker.io/onuroztornaci/methylhead-pipeline:qc_container`  | R 4.4.1, Quarto                 |

Build your own images → see [`/container-def-files`](/container-def-files/readme.md).

---

## 8 · Bundled panel and target files

* `data/blood-cell-type-reference.csv.gz` — Cell-type-specific reference for cell-count estimation
* `input/panel.csv` — Targeted CpG coordinates

Override with `--cell_reference` and `--panel` if you have a different panel.

---

## 9 · Utilities

### Extracting a basic dataset from [MethylKit](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-10-r87) output files

In cases where reads have already been aligned using another pipeline
and DNA methylation levels have been extracted using MethylKit,
we provide the script
[scripts/methylkit-to-matrix.sh](scripts/methylkit-to-matrix.sh)
to merge these outputs into an analyzable dataset
including the following output files:

* `cell-counts.csv` Cell count estimates assuming that the DNA methylation was measured in peripheral blood
* `coverage.csv` A matrix of read counts per sample per CpG site
* `illumina.csv` A matrix of methylation levels per sample per CpG site contained on the Illumina 450K or EPIC beadchip
* `methylation-scores.csv` A matrix of published methylation scores per sample calculated from the methylation levels (see [meffonym](https://github.com/perishky/meffonym) for details)
* `methylation.csv` A matrix of methylation levels per sample per CpG site for all sites contained in at least 10 reads for at least 50% of the samples
* `qc.html` A report providing an overview of the extracted data

The script can be executed as follows:

```
REPO=path/to/this/repository

$REPO/scripts/methylkit-to-matrix.sh \
    $REPO \
	path/to/data \
	path/to/target/panel \
	path/to/output \
	[hg19|hg38]
```

The final argument for the script, the genome assembly is optional.
The default is 'hg19'.  If read alignment was to the hg38 assembly, then
'hg38' needs to be specified for the DNA methylation score calculations
to be performed correctly.

This script assumes that R is installed
with the following R packages:

* IlluminaHumanMethylation450kanno.ilmn12.hg19
* IlluminaHumanMethylationEPICanno.ilm10b4.hg19
* IlluminaHumanMethylationEPICv2anno.20a1.hg38
* data.table
* dplyr
* methylKit
* [meffonym](https://github.com/perishky/meffonym)

---

## 10 · Troubleshooting cheatsheet

| Symptom                       | Likely cause & fix                                                        |
| ----------------------------- | ------------------------------------------------------------------------- |
| `ERROR: Apptainer not found`  | Install Apptainer ≥ 1.1 and add it to `$PATH`.                            |
| Java `<11` warning            | Forgot to `conda activate methylhead`.                                    |
| `No FASTQ files`              | Check your `--data` glob – must end in `.fastq.gz`.                       |
| `Index not found for hg19.fa` | Run **4 · reference build** or point `--genome_folder` to an indexed ref. |
| Path not mounted: data/reference outside `$HOME` | Move data and reference folders inside `$HOME`, **or** start Apptainer with `-B /abs/path:/abs/path` to bind-mount them. |
---

Happy methylating 🧬🚀
