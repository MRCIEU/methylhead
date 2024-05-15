## dnam-lung-cancer-pipeline

For nextflow usage Description:

Pipeline for analysing data generated from the DNAm lung cancer screening panel (https://github.com/MRCIEU/dnam-lung-cancer-screening-panel).

### Conda setup

Dependencies:

    bcftools     1.10         
    bedtools     2.31.1
    bismark      0.24.2        
    bowtie2      2.3.5.1
    bwa          0.7.18 
    bwameth      0.2.7
    fastqc       0.12.1        
    hisat2       2.2.0         
    multiqc      1.21 
    methyldackel 0.6.1
    nextflow     23.10.1
    perl         5.32.1 
    picard       2.18.23 
    qt           5.6.3         
    rstudio      1.1.456
    sambamba     1.0
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
wget https://www.twistbioscience.com/sites/default/files/resources/2022-06/covered_targets_Twist_Methylome_hg19_annotated_collapsed_final.bed.zip
mv covered_targets_Twist_Methylome_hg19_annotated_collapsed_final.bed.zip ${GENOMES}
samtools faidx hg19.fa
picard CreateSequenceDictionary REFERENCE=hg19.fa OUTPUT=hg19.fa.dict
picard BedToIntervalList \
I=covered_targets_Twist_Methylome_hg19_annotated_collapsed_final.bed \
O=covered_targets_Twist_Methylome_hg19_annotated_collapsed_final \
SD=hg19.fa.dict

```
Running time is about 2 hours.

## Prepare Beta Cell Types for Estimating cell counts

**The script reads blood panel data from a CSV file*

**It writes the unique blood panel data to a CSV file*

**The script iterates through each row of the blood panel to extract genomic regions*

**It uses the 'wgbstools' command-line tool to extract data for each region and appends it to a file*

**After that, it processes the extracted data to create a BED file containing relevant columns*

**Finally, it converts beta values to table format using 'wgbstools' and saves the result to a CSV file*

**We used the 'wgbstools' tool obtained from LoYfer for data extraction and conversion*

```
# This script prepares beta cell types files for analysis.

# Read blood panel data
blood_panel <- fread("panel-reduced_BS.csv")

# Remove rows where start equals end
blood_panel <- subset(blood_panel, blood_panel$start != blood_panel$end)

# Write unique blood panel data to a CSV file
write.csv(blood_panel, "blood_panel_unique.csv", row.names = FALSE)

# Assign blood panel to another variable
panel <- blood_panel

# Iterate through each row of the panel to extract regions
for (i in 1:nrow(panel)) {
  # Define region
  region <- paste0("chr", panel$chr[i], ":", panel$start[i], "-", panel$end[i])
  
  # Execute system command to extract data and append to file
  system(paste("./wgbstools convert -r", region, "| head -n1 >> panel-reduced-BS-sites.txt"))
}

# Extract relevant columns from the file and create a BED file
awk 'BEGIN {FS="[:, -]"; OFS="\t"} {start=$1+$2; end=$3+$4; print $1, start,end , $11,$12}' panel-reduced-BS-sites.txt | tail -n +1 > blood_cell.bed

# Convert beta values to table format and save to CSV file
./wgbstools beta_to_table blood_cell.bed --betas *.beta | column -t > blood_cell_types.csv
```


## Running the pipeline

**Bismark Usage:**

```
nextflow Bismark_pipeline.nf \
  --reads "[fastq path]/*_R{1,2}*.fastq.gz" \
  --t_param 4 \
  --memory_param 10000 \
  --genome_folder [genome index path]
  --multicore 4 \
  preview -with-dag flowchart.html \
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

**Picard Usage:**

```

nextflow Picard_pipeline.nf \
  --reads "[fastq path]/*_R{1,2}*.fastq.gz" \
  --intervals covered_targets_Twist_Methylome_hg19_annotated_collapsed_final
  --genome_folder  [BWA genome index path]

```


**Instructions:**

Input files (file1_1.fastq.gz, file1_2.fastq.gz, etc.), the main.nf script, and the output directory (outdir) should all be located in the same folder.

For parallel computing, organize the FASTQ files into separate folders. Each folder represents a distinct sequencing job.

Example 

folder1 [file1_1.fastq.gz, file1_2.fastq.gz, file2_1.fastq.gz, file2_2.fastq.gz]

folder2 [file3_1.fastq.gz, file3_2.fastq.gz, file3_1.fastq.gz, file3_2.fastq.gz]

folder3 [file4_1.fastq.gz, file4_2.fastq.gz, file5_1.fastq.gz, file5_2.fastq.gz]

this would be 3 sequencing job. 

All results will be output to the same directory specified by --outdir "Bismark_Results (or Picard_Results)"


## To do

* Genome indexing should be a process that runs only if the genome index needs to be created.
```

* Note: It is now possible to have a single directory with all fastq files.
  If new files are generated, just copy them to the directory and
  rerun the pipeline with the '-resume' option
