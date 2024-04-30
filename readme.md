## dnam-lung-cancer-pipeline

For nextflow usage Description:

Pipeline for analysing data generated from the DNAm lung cancer screening panel (https://github.com/MRCIEU/dnam-lung-cancer-screening-panel).

### Conda setup

Dependencies:

    bcftools     1.10         
    bedtools     2.31.1        
    bismark      0.24.2        
    bowtie2      2.3.5.1             
    fastqc       0.12.1        
    hisat2       2.2.0         
    multiqc      1.21          
    perl         5.32.1        
    qt           5.6.3         
    rstudio      1.1.456       
    samtools     1.18          
    trim-galore  0.6.10        
    trimmomatic  0.39          

To install these, you will need the following conda channels:
  - conda-forge
  - bioconda
  - defaults

Check which channels you have:
```
conda config --show channels
```

Add any missing channels like this:
```
conda config --add channels bioconda
```

Create a conda environment for the analysis and install packages
```
conda create -n Bismark python=3.8
pip install multiqc
conda install -c conda-forge trim-galore
conda install -c bioconda nextflow bismark samtools trimmomatic fastqc bedtools cutadapt bowtie
conda install -c conda-forge python-isal
```

## Prepare reference genome

```
GENOMES=[path to converted genome index]

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
mv hg19.fa.gz ${GENOMES}
bismark_genome_preparation --path_to_aligner /data/ms13525/lib/miniconda2/envs/Bismark/bin  --verbose ${GENOMES}
```

Running time is about 2 hours.

## Running the pipeline

**Usage:**

```
nextflow nextflow.nf \
  --reads "[fastq path]/*_R{1,2}*.fastq.gz" \
  --t_param 4 \
  --memory_param 10000 \
  --genome_folder [genome index path]
  --multicore 4 \
  preview -with-dag flowchart.html \
  --outdir "results" \
  -resume \
  -with-timeline time_line.html \
  -with-report report.html
```

**Parameters:**

- u_param (Bismark)

  This parameter determines how many reads will be used. 

- t_param (Fastqc)

  Fastqc Specifies the number of files which can be processed
  simultaneously. Each thread will be allocated 250MB of memory so you
  shouldn't run more threads than your available memory will cope
  with, and not more than 6 threads on a 32 bit machine.

- memory_param (Fastqc)

  Sets the base amount of memory, in Megabytes, used to process each
  file. Defaults to 512MB. You may need to increase this if you have a
  file with very long sequences in it. Allowed range (100 - 10000)

- multicore (Bismark)

  Number of cores to be used for Bismark Alignement. Allowed range (1-8)

- cores (Trim Galore!)

  Number of cores to be used for Trimming. It seems that --cores 4 could be a sweet spot, anything above has diminishing returns.    

 preview -with-dag flowchart
 
 The pipeline will be represented as a direct acyclic graph (DAG)

- with-timeline
 
 Using to enable the creation of the timeline report.

- with-report
 
 It creates an HTML execution report: a single document about resources usage (which includes many useful metrics about a workflow execution).


**Instructions:**

Input files (file1_1.fastq.gz, file1_2.fastq.gz, etc.), the main.nf script, and the output directory (outdir) should all be located in the same folder.

For parallel computing, organize the FASTQ files into separate folders. Each folder represents a distinct sequencing job.

Example 

folder1 [file1_1.fastq.gz, file1_2.fastq.gz, file2_1.fastq.gz, file2_2.fastq.gz]

folder2 [file3_1.fastq.gz, file3_2.fastq.gz, file3_1.fastq.gz, file3_2.fastq.gz]

folder3 [file4_1.fastq.gz, file4_2.fastq.gz, file5_1.fastq.gz, file5_2.fastq.gz]

this would be 3 sequencing job. 

All results will be output to the same directory specified by --outdir "results"


## To do

* Genome indexing should be a process that runs only if the genome index needs to be created

* Add processes for steps that follow, e.g.

```
cd results
bismark2report
bismark2summary
multiqc . -x .nextflow
Rscript R_codes.R
```

* 'nextflow.nf' should be renamed 'main.nf' (seems to be the convention for nextflow repos)

* 'R_codes.R' should be renamed 'generate-report.r' or something similar

* Note: It is now possible to have a single directory with all fastq files.
  If new files are generated, just copy them to the directory and
  rerun the pipeline with the '-resume' option

* 'workflow' code could be simplified

* for some reason the fastqc outputs got to 'results' and the
  remaining go to 'results/results', possibly the ${params.outdir} in
  the process outputs paths s unnecessary
