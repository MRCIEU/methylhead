This file, contains flowcharts to visualize and understand the workflow of both bismark and picard pipeline.

## Picard Pipeline Steps
* FastQC 

Quality control checks on raw sequence data to generate quality metrics.

Input  : fastq files

Output : fastq and html files

* Trim Galore

Trimming adapter sequences and low-quality bases from raw sequence reads.

Input  : fastq files

Output : fq.gz  files

* Alignment

Aligning trimmed sequence reads to a reference genome.

Input  : fq.gz files

Output : bam   files

* Sambamba 

Processing and manipulating aligned sequence data, such as sorting and indexing BAM files.

Input  : bam files

output : bam files 

* Sorted BAM Files

Generating sorted BAM files from the aligned data for downstream analysis.

Input  : bam files

output : bam files 

* Mark Duplicates 
Identifying and marking duplicate reads in the BAM files.
Input  : bam files
output : bam files and bam.bai files

* Collect HS Metrics 
Collecting hybrid selection metrics, which measure the efficiency of target enrichment.
Input  : bam files
Output : txt files and covarege files

* Collect MM Metrics 
Collecting mismatch metrics to assess alignment accuracy.
Input  : bam files
Output : alignment_summary_metrics, gc_bias_detail_metrics , gc_bias(.pdf), summary_metrics, insert_size_histogram.pdf , insert_size_metrics

* MethylDackel 
Performing methylation calling to determine the methylation status of cytosines.
Input  : bam files
Output : svg files  

* bedGraph
Generating bedGraph files, which are text files describing genome-wide data in a graph format.
Input  : bam files files
Output : bedGraph files

* Processed bedGraph
Further processing of bedGraph files for downstream analysis.
Input  : bedGraph files
Output : bedGraph files

* Samtools Stats
Generating statistics from BAM files using Samtools.
Input  : bam files
Output : samtools_stats(.txt), markdup_samtools_stats(.txt) files

* MethylKit 
Analyzing differential methylation using the MethylKit tool.
Input  : bam files
Output : methylKit files

* Methylation Matrix
Creating a matrix of methylation values for each sample.
Input  : methylKit files
Output : csv and pdf files

* DNAm Full Matrix
Generating a comprehensive DNA methylation matrix.
Input  : methylKit files
Output : csv files

* Estimate Cell Counts
Estimating cell type composition based on methylation data.
Input  : csv files
Output : csv files

* DNA Methylation Scores 
Calculating scores that reflect the methylation status of specific regions.
Input  : csv files
Output : csv files

* MultiQC
Aggregating and visualizing quality control metrics from multiple sources.
Input  : txt files
Output : html files

## Bismark Pipeline Steps
* FastQC
Quality control checks on raw sequence data to generate quality metrics.
Input  : fastqc files
Output : fastqc files

* Trim Galore
Trimming adapter sequences and low-quality bases from raw sequence reads.
Input  : fastqc files
Output : fq.gz files

* Alignment
Aligning trimmed sequence reads to a reference genome.
Input  : fq.gz files
Output : bam files

* Deduplication
Identifying and marking duplicate reads in the BAM files.
Input  : bam files
Output : bam files files

* Methylation Extraction
Performing methylation calling to determine the methylation status of cytosines.
Input  : bam files
Output : deduplicated.bedGraph.gz, bismark.cov.gz, deduplicated_splitting_report.txt, deduplicated.M-bias.txt

* Reports
Compiling various reports generated throughout the pipeline.
Input  : Various intermediate files
Output : Report files (various formats)

* Methylation Matrix
Creating a matrix of methylation values for each sample.
Input  : methylKit files
Output : csv and pdf files

* DNAm Full Matrix
Generating a comprehensive DNA methylation matrix.
Input  : methylKit files
Output : csv files

* Estimate Cell Counts
Estimating cell type composition based on methylation data.
Input  : csv files
Output : csv files

* DNA Methylation Scores 
Calculating scores that reflect the methylation status of specific regions.
Input  : csv files
Output : csv files

* MultiQC
Aggregating and visualizing quality control metrics from multiple sources.
Input  : txt files
Output : html files
