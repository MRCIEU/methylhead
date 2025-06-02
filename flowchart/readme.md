## Overview

- This workflow processes Whole Genome Bisulfite Sequencing (WGBS) raw data through a series of bioinformatics steps for DNA methylation analysis, quality control, and association testing. 

- It integrates alignment, methylation calling, QC metrics extraction, and advanced scoring methods to generate statistically interpretable methylation data matrices.

---

## Detailed Steps with Inputs, Outputs, and Data Characteristics

### 1. Raw FASTQ Reads  
- **Input Data:** Paired-end or single-end raw sequencing reads in FASTQ format (`*.fastq.gz`).  
- **Purpose:** Provide raw sequence data for downstream trimming, alignment, and quality control.

### 2. Trim Galore  
- **Input:** Raw FASTQ reads  
- **Output:** Trimmed FASTQ reads with adapters and low-quality bases removed.  
- **Method:** Adapter trimming + quality filtering (Phred score cutoff typically 20-30).  
- **Effect:** Improves read quality and alignment accuracy by removing sequencing artifacts.

### 3. Alignment - Bwa-meth  
- **Input:** Trimmed FASTQ reads  
- **Output:** BAM files with reads aligned to the reference genome using BWA-meth (bisulfite-aware aligner).  
- **Characteristics:** Produces sorted, indexed BAM files with methylation-aware alignment flags.  
- **Statistical Note:** Alignment rate and mapping quality scores are essential QC metrics for downstream filtering.

### 4. Sambamba Sort  
- **Input:** Raw BAM from alignment  
- **Output:** Position-sorted BAM file  
- **Purpose:** Necessary for efficient duplicate marking and indexing.

### 5. Mark Duplicates  
- **Input:** Sorted BAM  
- **Output:** BAM with PCR duplicates flagged  
- **Statistical Rationale:** Duplicate reads can bias methylation estimates; marking duplicates prevents overcounting.

---

### 6. MethylDackel  
- **Input:** Duplicate-marked BAM  
- **Output:** Per-CpG methylation calls, including methylated/unmethylated read counts, methylation percentage per site.  
- **Data Format:** TSV with columns for chromosome, position, methylated reads, unmethylated reads, and methylation percentage.  
- **Statistical Use:** Basis for estimating methylation beta values and confidence intervals.

### 7. BedGraph Generation  
- **Input:** Duplicate-marked BAM  
- **Output:** BedGraph files representing genome-wide methylation coverage and levels.  
- **Usage:** Visual inspection of methylation landscapes in genome browsers.

### 8. Processed BedGraph  
- **Input:** Raw BedGraph  
- **Output:** Filtered and normalized BedGraph files  
- **Purpose:** Remove low coverage or low-quality sites to improve downstream analyses.

### 9. Samtools Stats  
- **Input:** BAM files  
- **Output:** Summary statistics (e.g., total reads, mapped reads, duplicates)  
- **Use:** Monitor data quality and sequencing depth.

---

### 10. Collect HS Metrics  
- **Input:** BAM files  
- **Output:** Hybrid selection (target enrichment) specific QC metrics such as on-target rate, fold enrichment.  
- **Relevance:** Validates target capture efficiency in targeted bisulfite sequencing.

### 11. Collect MM Metrics  
- **Input:** BAM files  
- **Output:** Methylation-specific QC metrics like bisulfite conversion rate.  
- **Importance:** Ensures experimental quality of methylation data.

---

### 12. MethylKit Analysis  
- **Input:** MethylDackel output files (per-sample methylation calls)  
- **Output:** Aggregated methylation matrix (CpG sites × samples) with methylation percentages (beta values).  
- **Data Format:** Matrix with CpG site identifiers (chromosome_position) as rows and samples as columns.  
- **Statistical Note:** Beta values range from 0 to 1, representing methylation proportion.  
- **Use:** Basis for differential methylation analysis and epigenome-wide association studies (EWAS).

### 13. DNA Methylation Scores  
- **Input:** Methylation matrix  
- **Output:** Composite methylation scores per sample or genomic region (e.g., average beta values or principal components).  
- **Statistical Role:** Summarizes methylation variability, reduces dimensionality for modeling.

### 14. Estimate Cell Counts  
- **Input:** Methylation matrix  
- **Output:** Estimated proportions of cell types (e.g., blood cell fractions) using reference-based deconvolution methods.  
- **Statistical Model:** Regression or constrained projection methods applied to methylation signatures.

### 15. Illumina Matrix  
- **Input:** Methylation matrix (platform-specific)  
- **Output:** Normalized matrix tailored for Illumina arrays, including batch-corrected beta values.  
- **Purpose:** Compatibility with downstream statistical tools for Illumina data.

### 16. Coverage Matrix  
- **Input:** BAM files  
- **Output:** Matrix of read depths per CpG site per sample.  
- **Use:** Quality filtering thresholding (e.g., removing sites with coverage <10).

---

### 17. BSmap Alignment  
- **Input:** Raw FASTQ reads  
- **Output:** BAM files aligned with BSmap bisulfite aligner (alternative method).  
- **Comparison:** Used to validate alignment results and downstream methylation calls.

### 18. CAMDA Calculation  
- **Input:** BSmap BAM files  
- **Output:** CAMDA concurrence scores measuring simultaneous methylation/demethylation events.  
- **Statistical Purpose:** Identify CpG sites with coordinated epigenetic changes.

### 19. CAMDA Output  
- **Input:** CAMDA calculation results  
- **Output:** Processed CAMDA scores ready for integration.

### 20. CAMDA Matrix  
- **Input:** CAMDA output  
- **Output:** Matrix of CAMDA scores (CpG sites × samples) for association testing.

---

### 21. MultiQC  
- **Input:** QC outputs: FastQC reports, HS Metrics, MM Metrics, Samtools Stats, Illumina Matrix, Cell Counts, Methylation Scores, Coverage Matrix, CAMDA Matrix  
- **Output:** Comprehensive multi-tool QC report (HTML).  
- **Use:** Centralized overview of sequencing and methylation data quality.

### 22. QC Report  
- **Input:** MultiQC summary  
- **Output:** Final QC assessment document.

### 23. Association Test  
- **Input:** QC-passed methylation data and phenotype data  
- **Output:** Statistical tests (e.g., linear regression, mixed models) results associating methylation levels with phenotypes.  
- **Statistical Methods:** Adjustment for confounders (age, sex, cell counts), multiple testing correction (FDR).  
- **Result Format:** Tables of CpG sites with p-values, effect sizes, confidence intervals.

---
