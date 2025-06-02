# Minimal Test Dataset

- This directory provides a **small, reproducible test dataset** derived from real WGBS data for fast pipeline prototyping and demonstration. 
- All files are based on a subset of the original study to keep size and runtime minimal.

## Content Overview

- **test-data/**: Mini paired-end FASTQ files (subsampled from ENA [PRJNA730913](https://www.ebi.ac.uk/ena/browser/view/PRJNA730913), filtered for selected regions).
- **test-reference/**: Mini hg19 reference FASTA, containing only test regions.
- **test.bed**: BED file for target regions.
- **test-panel.csv**: Example design file listing loci and regions.
- **test-phenotype.csv**: Example phenotype data for test samples.
- **test-models.csv**: Example model parameters.
- **bash-files/**: All-in-one scripts to (re)generate test data and reference from public resources.

> All files can be **fully regenerated** using the provided scripts.  
> No private or protected data is included – everything is based on public and open-access resources.

