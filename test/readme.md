## Minimal Test Dataset

- Includes a **small, reproducible dataset** derived from real WGBS data ([ENA PRJNA730913](https://www.ebi.ac.uk/ena/browser/view/PRJNA730913)), designed for fast prototyping and demonstration.
- All files are based on a subset of the original study, selected to minimize runtime and disk usage.
- Every file can be **fully regenerated** using the provided scripts.  
- No private or protected data is included — everything is public and open-access.

---
### Content Overview

- **test-data/**: Mini paired-end FASTQ files (subsampled and region-filtered from ENA).
- **test-reference/**: Minimal hg19 reference FASTA, containing only test regions.
- **test-target.bed**: BED file with target regions.
- **test-panel.csv**: Example panel design file listing loci and regions.
- **test-phenotype.csv**: Example phenotype data for test samples.
- **test-models.csv**: Example model parameters for the pipeline.
- **bash-files/**: All-in-one scripts to (re)generate test data, references, and matrices from public sources.

---
## Quick Start

Follow these steps to generate all required files for a minimal methylation analysis demo:

1. **Download FASTQ files**  
   Download public test FASTQ files:  
   ```
   bash bash-files/fastq-files-download.sh
   ```

2. **Align FASTQ files with BWAmeth and extract region-specific paired-end FASTQ files**  
   - Align all downloaded FASTQ files to the reference genome (`hg19.fa`) using BWAmeth.  
   - Generate small, region-specific (blood-cell-types-regions.bed) paired-end FASTQ files in the `test-data` directory:  
   ```
   bash bash-files/create-example-data.sh <bed-file> <output-dir>
   ```

3. **Generate extended BED regions and subset FASTQs**  
   - Using the indexed `hg19.fa` and the new small FASTQ files as input, run the pipeline to generate:  
     - Sorted BAM files  
     - CpG-level methylation matrix  
     - Illumina-450k beta-value matrix  
     - The main output, `methylation-illumina-regions.bed` (the union of all regions covered in the two matrices), is provided in this repository.  
   ```
   bash bash-files/create-regions.sh <bam-files/*.bam> <*.bed>
   ```

4. **Create minimal reference FASTA**  
   - Create a `test-reference.fa`, a minimal reference FASTA file for the selected regions, and index it with BWAmeth.  
   - Create a `test-target.bed` file.  
   ```
   bash bash-files/reference-create.sh
   ```
---
## Scripts in bash-files/

- `fastq-files-download.sh`  
  Downloads public test FASTQ files (step 1).
- `create-example-data.sh`  
  Aligns all FASTQ files to the reference genome (`hg19.fa`) with BWAmeth and generates small, region-specific paired-end FASTQ files (step 2).
- `create-regions.sh`  
  Using the indexed `hg19.fa` and the new small FASTQ files as input, generates sorted BAM files, CpG-level methylation matrix, Illumina-450k beta-value matrix, and outputs `methylation-illumina-regions.bed` (step 3).
- `reference-create.sh`  
  Creates a minimal reference FASTA (`test-reference.fa`) and its index, and also generates `test-target.bed` for the selected regions (step 4).

---
## Requirements

- **BWAmeth**
- **BEDtools**
- **Samtools**
- **Python 3**
- **hg19 reference genome** (indexed for BWAmeth and Samtools)
- Standard UNIX tools: `bash`, `awk`
---
## Notes

- All steps are fully reproducible and use only public data.
- Adjust script paths or filenames as needed for your own environment.
- No private or restricted data is used.
---
