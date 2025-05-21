# Workflow Overview

This file contains flowcharts to visualize and understand the workflow of both the Bismark and Picard pipelines.

## Picard Pipeline Steps

### FastQC
**Description**: Quality control checks on raw sequence data to generate quality metrics.  
**Input**: fastq files  
**Output**: fastq and html files  

### Trim Galore
**Description**: Trimming adapter sequences and low-quality bases from raw sequence reads.  
**Input**: fastq files  
**Output**: fq.gz files  

### Alignment
**Description**: Aligning trimmed sequence reads to a reference genome.  
**Input**: fq.gz files  
**Output**: bam files  

### Sambamba
**Description**: Processing and manipulating aligned sequence data, such as sorting and indexing BAM files.  
**Input**: bam files  
**Output**: bam files  

### Sorted BAM Files
**Description**: Generating sorted BAM files from the aligned data for downstream analysis.  
**Input**: bam files  
**Output**: bam files  

### Mark Duplicates
**Description**: Identifying and marking duplicate reads in the BAM files.  
**Input**: bam files  
**Output**: bam files and bam.bai files  

### Collect HS Metrics
**Description**: Collecting hybrid selection metrics, which measure the efficiency of target enrichment.  
**Input**: bam files  
**Output**: txt files and coverage files  

### Interval Files
**Description**: Interval files from the panel for collecting HS Matrix.  
**Input**: bed file  
**Output**: interval file  

### Collect MM Metrics
**Description**: Collecting mismatch metrics to assess alignment accuracy.  
**Input**: bam files  
**Output**: alignment_summary_metrics, gc_bias_detail_metrics, gc_bias (.pdf), summary_metrics, insert_size_histogram.pdf, insert_size_metrics  

### MethylDackel
**Description**: Performing methylation calling to determine the methylation status of cytosines.  
**Input**: bam files  
**Output**: svg files  

### bedGraph
**Description**: Generating bedGraph files, which are text files describing genome-wide data in a graph format.  
**Input**: bam files  
**Output**: bedGraph files  

### Processed bedGraph
**Description**: Further processing of bedGraph files for downstream analysis.  
**Input**: bedGraph files  
**Output**: bedGraph files  

### Samtools Stats
**Description**: Generating statistics from BAM files using Samtools.  
**Input**: bam files  
**Output**: samtools_stats (.txt), markdup_samtools_stats (.txt) files  

### MethylKit
**Description**: Analyzing differential methylation using the MethylKit tool.  
**Input**: bam files  
**Output**: methylKit files  

### Illumina (Methylation) Matrix
**Description**: Creating a matrix of methylation values for each sample and CpGs sites. 
**Input**: methylKit files  
**Output**: csv and pdf files  

### DNAm Matrix
**Description**: Generating a comprehensive DNA methylation matrix.  
**Input**: methylKit files  
**Output**: csv files  

### Estimate Cell Counts
**Description**: Estimating cell type composition based on methylation data.  
**Input**: csv files  
**Output**: csv files  

### DNA Methylation Scores
**Description**: Calculating scores that reflect the methylation status of specific regions.  
**Input**: csv files  
**Output**: csv files  

### MultiQC
**Description**: Aggregating and visualizing quality control metrics from multiple sources.  
**Input**: txt files  
**Output**: html files  

