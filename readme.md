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
conda install -c bioconda nextflow bismark samtools trimmomatic fastqc bedtools cutadapt bowtie bwa bwameth sambamba methyldackel
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

## Usage: 

```
nextflow main.nf --pipeline bismark \
--reads "[fastq path]/*_R{1,2}*.fastq.gz" \
--genome_folder [genome index path] \
-resume

nextflow main.nf --pipeline picard \
--reads "[fastq path]/*_R{1,2}*.fastq.gz" \
--intervals covered_targets_Twist_Methylome_hg19_annotated_collapsed_final \
--genome_folder [BWA genome index path] \
-resume 
```

*Camda Pipeline*

```
nextflow main.nf --reads "[fastq path]/*_R{1,2}*.fastq.gz" \ 
--genome_folder [genome folder path] \
--samtools_path [samtools fodler path] \
-resume
```

**Instructions:**

Input files (file1_1.fastq.gz, file1_2.fastq.gz, etc.), the main.nf script, and the output directory (outdir) should all be located in the same folder.

For parallel computing, organize the FASTQ files into separate folders. Each folder represents a distinct sequencing job.

Example 

folder1 [file1_1.fastq.gz, file1_2.fastq.gz, file2_1.fastq.gz, file2_2.fastq.gz]

folder2 [file3_1.fastq.gz, file3_2.fastq.gz, file3_1.fastq.gz, file3_2.fastq.gz]

folder3 [file4_1.fastq.gz, file4_2.fastq.gz, file5_1.fastq.gz, file5_2.fastq.gz]

this would be 3 sequencing job. 


## To do

* Genome indexing should be a process that runs only if the genome index needs to be created.
```

* Note: It is now possible to have a single directory with all fastq files.
  If new files are generated, just copy them to the directory and
  rerun the pipeline with the '-resume' option
