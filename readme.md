# **Methylhead** – DNAm Lung‑Cancer‑Screening Pipeline

**Methylhead** is a modular **Nextflow** workflow that turns raw FASTQ files from the DNAm lung‑cancer‑screening panel into QC‑checked methylation matrices, cell‑composition estimates, and model‑based risk scores—ready for downstream statistics or reporting.

> **Default reference files** for the panel are already shipped with the repo:
>
> * `data/blood_cell_types_extended.bed` – CpG coordinates extended with blood‑cell reference loci used for **cell‑composition estimation**
> * `input/panel.csv` – per‑locus metadata & QC thresholds
>
> Supply your own via `--panel` / `--panel_qc` if you work with a different panel.

Panel manifest and documentation: [https://github.com/MRCIEU/dnam-lung-cancer-screening-panel](https://github.com/MRCIEU/dnam-lung-cancer-screening-panel)

---

## 1 · Install

```bash
# Conda env with Nextflow & Java ≥ 11
conda create -n methylhead nextflow -c bioconda
conda activate methylhead
```

---

## 2 · Prepare reference genome (one‑off)

```bash
bash prepare-reference-genome.sh -N you@example.com  # ▸ writes to reference/ 
```
>-N Address for auto‑emailed report (optional)

---

## 3 · Quick start (public FASTQ + simulated phenotypes)

Follow the steps below to fetch the demo data, execute the workflow, and review the results.

1. **Download data** – 20 real paired‑end FASTQ files from ENA study \[PRJNA730913] are placed in `test-data/`.

   ```bash
   bash test-data.sh
   ```
2. **Run the DNAm‑panel workflow**

   ```bash
   nextflow run main.nf -C nextflow-test.config -N you@example.com --resume
   ```

---

## 4 · Run on your own samples

```bash
nextflow run main.nf \
    --data            path/to/fastqs/ \
    --genome_folder   path/to/hg19.fa \
    --panel           path/to/panel.bed \
    --panel_qc        path/to/panel_qc.csv \
    --phenotype       path/to/phenotype.csv \
    --models          path/to/models.csv \
    --outdir          results/ \
    -N you@example.com \
    --resume
```

All CLI parameters can instead be written into a config and supplied with `-C nextflow.config`.

> **Note:** A demo models file (`input/models-test.csv`) lives in the `input/` folder. Edit this file—or point `--models` to your own CSV in `input/`—to run custom EWAS or risk‑prediction models.

### Key parameters

| Param           | Purpose                                    | Default (demo)                       |
| --------------- | ------------------------------------------ | ------------------------------------ |
| `data`          | Folder (or glob) of paired‑end gz FASTQs   | `test-data/`                         |
| `genome_folder` | Indexed reference FASTA                    | `reference/hg19.fa`                  |
| `panel`         | BED of CpG loci                            | `data/blood_cell_types_extended.bed` |
| `panel_qc`      | CSV with per‑locus thresholds              | `input/panel.csv`                    |
| `phenotype`     | Sample metadata (CSV)                      | `input/phenotype-test.csv`           |
| `models`        | EWAS / risk‑prediction models (CSV)        | `input/models-test.csv`              |
| `outdir`        | Where to write results                     | `results-test/`                      |
| `email` / `-N`  | Address for auto‑emailed report (optional) | *(none)*                             |

---

## 5 · Outputs

```
results/
├── alignments/          # dedup BAM + stats
├── methylation_calls/   # per‑sample BedGraphs & bigWigs
├── matrices/            # CpG, coverage & 450k matrices (TSV)
├── qc/                  # MultiQC + html/pdf report
└── predictions/         # Risk scores & association test results
```

---

## 6 · Reproducing the DAG

A pre‑generated pipeline graph [`workflow.png`](https://github.com/MRCIEU/dnam-lung-cancer-pipeline/blob/main/flowcharts/workflow.png) is committed to the `flowcharts/` directory, so you can inspect the workflow without running anything.

```bash
nextflow run main.nf -C … --resume -with-dag flow.svg
```

## 7 · Container images (automatic)

This workflow is shipped with three pre-built OCI/Apptainer images.  
When you launch the pipeline **Nextflow pulls them on-demand** (via
`oras://`) and attaches the right image to each process, so you don’t
have to install any tool chain manually.

| Logical image | Default URI                                                                                 | What it contains                              |
| ------------- | ------------------------------------------------------------------------------------------- | --------------------------------------------- |
| `wgbs_image`  | `oras://docker.io/onuroztornaci/methylhead-pipeline:wgbs_analysis`                          | WGBS aligners + core QC tools                 |
| `meth_image`  | `oras://docker.io/onuroztornaci/methylhead-pipeline:meth_analysis`                          | R 4.4.3 + Python + Bioconductor methylation   |
| `qc_image`    | `oras://docker.io/onuroztornaci/methylhead-pipeline:qc_container`                           | R 4.4.1 + Quarto for report generation        |

Override them on the CLI (e.g. `--wgbs_image my.registry/wgbs:1.2`) or in
a config file.

> **Prefer to build your own images?**  
> The repo also ships Apptainer definitions that reproduce the three
> images one-to-one.  
> See [`container-def-files`](https://github.com/MRCIEU/dnam-lung-cancer-pipeline/tree/main/container-def-files) for the full recipes.
