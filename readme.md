# **Methylhead** – DNAm Lung‑Cancer‑Screening Pipeline


**Methylhead** is a modular **Nextflow** workflow that turns raw FASTQ files from the DNAm lung‑cancer‑screening panel into QC‑checked methylation matrices, cell‑composition estimates and model‑based risk scores – ready for downstream statistics or reporting.

> **Default reference files** for the wet‑lab panel are already shipped with the repo:
>
> * `data/blood_cell_types_extended.bed` – CpG coordinates extended with blood‑cell reference loci used for **cell‑composition estimation**
> * `input/panel.csv` – per‑locus metadata & QC thresholds
>
> Supply your own via `--panel` / `--panel_qc` if you work with a different panel.

Panel manifest and documentation: [https://github.com/MRCIEU/dnam-lung-cancer-screening-panel](https://github.com/MRCIEU/dnam-lung-cancer-screening-panel)

---

## 1 · Install

```bash
# 1. Conda env with Nextflow & Java ≥ 11
conda create -n methylhead nextflow -c bioconda
conda activate methylhead

# 2. (one‑off) grab and index hg19
bash prepare-reference-genome.sh   # ▸ writes to reference/
```

---

## 2 · Quick start (public data)

Follow the steps below to fetch the demo data, execute the workflow, and review the results:

```bash
# 1) Download the public demonstration dataset (20 paired‑end FASTQs → test-data/)
bash test-data.sh  
#   ↳ retrieves samples from ENA study [PRJNA730913](https://www.ebi.ac.uk/ena/browser/view/PRJNA730913) 
# 2) Run the DNAm‑panel workflow on the demo data
nextflow run main.nf -C nextflow-test.config -N you@example.com --resume
```


## 3 · Run on your own samples

```bash
nextflow run main.nf \
    --data            path/to/fastqs/ \
    --genome_folder   reference/hg19.fa \
    --panel           input/my_panel.bed \
    --phenotype       input/phenotype.csv \
    --models          input/models.csv \
    --panel_qc        input/panel_qc.csv \
    --outdir          results/ \
    -N you@example.com \
    --resume
```

All CLI parameters can instead be written into a config and supplied with `-C myrun.config`.

### Key parameters

| param           | purpose                                    | default (demo)                       |
| --------------- | ------------------------------------------ | ------------------------------------ |
| `data`          | folder (or glob) of paired‑end gz FASTQs   | `test-data/`                         |
| `genome_folder` | indexed reference FASTA                    | `reference/hg19.fa`                  |
| `panel`         | BED of CpG loci                            | `data/blood_cell_types_extended.bed` |
| `panel_qc`      | CSV with per‑locus thresholds              | `input/panel.csv`                    |
| `phenotype`     | sample metadata (CSV)                      | `input/phenotype-test.csv`           |
| `models`        | EWAS / risk‑prediction models (CSV)        | `input/models-test.csv`              |
| `outdir`        | where to write results                     | `results-test/`                      |
| `email/-N`      | address for auto‑emailed report (optional) | *(none)*                             |

---

## 4 · Outputs

```
results/
├── alignments/          # dedup BAM + stats
├── methylation_calls/   # per‑sample BedGraphs & bigWigs
├── matrices/            # CpG, coverage & 450k matrices (TSV)
├── qc/                  # MultiQC + html/pdf report
└── predictions/         # risk scores & association test results
```

---

## 5 · Reproducing the DAG

A pre‑generated pipeline graph `workflow.png` is already committed to the flowcharts/ directory of this repository, so you can inspect the workflow without running anything.

```bash
nextflow run main.nf -C … -resume -with-dag flow.svg
```

