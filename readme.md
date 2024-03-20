## dnam-lung-cancer-pipeline
For nextflow usage Description:

Pipeline for analysing data generated from the DNAm lung cancer screening panel (https://github.com/MRCIEU/dnam-lung-cancer-screening-panel).
conda create -n Bismark
conda install -c bioconda multiqc bismark trim-galore samtools trimmomatic fastqc bedtools cutadapt bs-seeker2 bowtie

name: Bismark
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:

    bcftools     1.10         
    bedtools     2.31.1        
    bismark      0.24.2        
    bowtie2      2.3.5.1       
    bs-seeker2   2.1.7         
    fastqc       0.12.1        
    hisat2       2.2.0         
    multiqc      1.8           
    perl         5.32.1        
    qt           5.6.3         
    rstudio      1.1.456       
    samtools     1.18          
    trim-galore  0.6.10        
    trimmomatic  0.39          

Refernce Genome       : hg38
u_param (Bismark)     : This parameter determines how many reads will be used.  
t_param (Fastqc)      : Fastqc Specifies the number of files which can be processed
                        simultaneously. Each thread will be allocated 250MB of
                        memory so you shouldn't run more threads than your
                        available memory will cope with, and not more than
                        6 threads on a 32 bit machine. 
memory_param (Fastqc) : Sets the base amount of memory, in Megabytes,  used to process
                        each file. Defaults to 512MB. You may need to increase this if
                        you have a file with very long sequences in it.
                        Allowed range (100 - 10000)

nextflow nextflow.nf --reads 'path/*_{1,2}.fastq.gz' --t_param number --u_param number --memory_param number --genome_folder 'bismark.ref path'

