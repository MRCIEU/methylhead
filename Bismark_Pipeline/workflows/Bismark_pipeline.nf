#!/usr/bin/env nextflow

log.info"""\
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  Nextflow DNA Methylation Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *** Steps          
     1. QC Steps
     2. Bismark Alignment& Methylation steps
     3. Methylation Matrix  
     4. Estimation cell counts
     5. DNA Methylation Score
           
 *** Parameters Info   
     Refernce Genome       : hg19 (or hg38)
     u_param (Bismark)     : This parameter determines how many reads will be used.  
     t_param (Fastqc)      : Fastqc Specifies the number of files which can be processed
                             simultaneously. Each thread will be allocated 250MB of
                             memory so you shouldn't run more threads than your
                             available memory will cope with, and not more than
                             6 threads on a 32 bit machine. 
     memory_param (Fastqc) : Sets the base amount of memory, in Megabytes,  used to process
                             each file. Defaults to 512MB. You may need to increase this if
                             you have a file with very long sequences in it. Allowed range (100 - 10000)
     multicore (Bismark)   : Number of cores to be used for Bismark Alignement. Allowed range (1-8)
     cores (Trim Galore!)  : Number of cores to be used for Trimming.
                             It seems that --cores 4 could be a sweet spot, anything above has diminishing returns.       
"""
include { Fastqc } from '../modules/Fastqc'
include { Trim_galore } from '../modules/Trim_galore'
include { Bismark_alignment } from '../modules/Bismark_alignment'
include { Bismark_Deduplication } from '../modules/Bismark_Deduplication'
include { Bismark_methylation_extraction } from '../modules/Bismark_methylation_extraction'
include { Methylation_Matrix } from '../modules/Methylation_Matrix'
include { Estimate_cell_counts } from '../modules/Estimate_cell_counts'
include { Bismark_scores } from '../modules/Bismark_scores'

workflow Bismark_pipeline {

    take:
    reads   
    outdir 
    
    main:
   read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)    
         t_param=params.t_param
         memory_param=params.memory_param
   Fastqc(read_pairs_ch,params.t_param,params.memory_param)
         cores=params.cores
   Trim_galore(read_pairs_ch,cores)
         trim_ch=Trim_galore.out.fq         
         u_param= params.u_param 
         multicore= params.multicore 
   bam_files_ch = Bismark_alignment(trim_ch,params.u_param,multicore)
        ch_bam = Bismark_alignment.out.bam
   Bismark_Deduplication(ch_bam) 
       dedup_bam=Bismark_Deduplication.out.bam
   Bismark_methylation_extraction(dedup_bam)
   coverage= Bismark_methylation_extraction.out.coverage               
        files_ch = coverage.collectFile(name:"*.cov.gz", newLine: true)
   Methylation_Matrix(files_ch)
   Estimate_cell_counts(files_ch)    
   Meth_Matrix = Methylation_Matrix.out
        files_ch2 = Meth_Matrix.collectFile(name:"*.csv", newLine: true)
   Bismark_scores(files_ch2)
}


log.info("""\
  +--------------------------------------------------------------------------------+
  | Pipeline Step               |         Description                              |
  +----------------------+---------------------------------------------------------+
  | Fastqc                      | Quality control of Fastq files                   |
  | Trimming                    | Removal of low-quality bases                     | 
  | Alignment                   | Alignment of paired-end fastq files to bam files |
  | Deduplication               | Removal of PCR duplicates                        |
  | DNA Methylation Extraction  | Analysis of DNA methylation                      |
  | Methylation Matrix          | Generation of methylation matrix                 |
  | Bismark Scores              | Calculation of DNA methylation indices           |
  +--------------------------------------------------------------------------------+
""")

workflow.onComplete {
    log.info("""\
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            THE END
                       Have a Great Day!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
""")
}
