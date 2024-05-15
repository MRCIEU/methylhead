#!/usr/bin/env nextflow

log.info"""\
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  Nextflow DNA Methylation Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *** Steps          
     1. QC Steps
     2. Bismark Alignment& Methylation steps
     3. Methylation Matrix 
     4. DNAm Full Matrix 
     5. Estimation cell counts
     6. DNA Methylation Score
     7. Multiqc 
           
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
include { Alignment } from '../modules/Alignment'
include { Deduplication } from '../modules/Deduplication'
include { Methylation_extraction } from '../modules/Methylation_extraction'
include { DNAm_Full_Matrix } from '../modules/DNAm_Full_Matrix'
include { Methylation_Matrix } from '../modules/Methylation_Matrix'
include { Estimate_cell_counts } from '../modules/Estimate_cell_counts'
include { DNA_Methylation_Scores } from '../modules/DNA_Methylation_Scores'
include { Reports } from '../modules/Reports'
include { Multiqc } from '../modules/Multiqc'

workflow Bismark_pipeline {

    take:
    reads   
    outdir 
    
    main:
     Fatqc_Files = Channel.fromFilePairs(params.reads, checkIfExists: true)    
         t_param=params.t_param
         memory_param=params.memory_param
     Fastqc(Fatqc_Files,t_param,memory_param)
         cores=params.cores
     Trim_galore(Fatqc_Files,cores)
         trim_ch=Trim_galore.out.fq         
         u_param= params.u_param 
         multicore= params.multicore 
     bam_files_ch = Alignment(trim_ch,u_param,multicore)
        ch_bam = Alignment.out.alignment_bam
     Deduplication(ch_bam) 
         dedup_bam=Deduplication.out.dedup_bam
     Methylation_extraction(dedup_bam)
     coverage= Methylation_extraction.out.coverage               
        files_ch = coverage.collectFile(name:"*.cov.gz", newLine: true)
     Methylation_Matrix(files_ch)
     DNAm_Full_Matrix(files_ch)
         full_matrix=DNAm_Full_Matrix.out
            full_matrix2=full_matrix.collectFile(name:"*.csv", newLine: true)
     Estimate_cell_counts(full_matrix2)    
     Meth_Matrix = Methylation_Matrix.out.meth_matrix
        files_ch2 = Meth_Matrix.collectFile(name:"*.csv", newLine: true)
     DNA_Methylation_Scores(files_ch2)
          alignment_text = Alignment.out.alignment_report
          deduplication_text = Deduplication.out.dedup_report
          methylation_text = Methylation_extraction.out.splitting_report 
          methylation_text2= Methylation_extraction.out.mbias
            alignment_t = alignment_text.collectFile(name:"*.txt",newLine:true)
            deduplication_t = deduplication_text.collectFile(name:"*.txt",newLine:true)
            methylation_t = methylation_text.collectFile(name:"*.txt",newLine:true)
            methylation_t2 = methylation_text2.collectFile(name:"*.txt",newLine:true)
     Reports(alignment_t,deduplication_t,methylation_t,methylation_t2) 
           Channel.empty()
          .mix( Fastqc.out )
          .mix( Trim_galore.out )
          .mix( Alignment.out )
          .mix( Deduplication.out )
          .mix( Methylation_extraction.out ) 
          .map { sample_id, files -> files }
          .collect()
          .set { log_files }  
     Multiqc(log_files)

}

log.info("""\
+---------------------------+----------------------------------------------=---+
| Pipeline Step             | Description                                      |
+---------------------------+--------------------------------------------------+
| QC Steps                  | Quality control of Fastq files                   |
| Trimming Fastqc Files     | Trimming and quality assessment of Fastq files   |
| Bismark Alignment         | Alignment of paired-end fastq files to bam files |
| Bismark Deduplication     | Removal of duplicate reads in Bismark alignment  |
| Methylation steps         | Creation of DNA methylation bedGraph             |
| Methylation Matrix        | Generation of methylation matrix (CpGs Sites)    |
| DNAm Full Matrix          | DNA methylation matrix                           |
| Estimation cell counts    | Estimation of cell counts from methylation data  | 
| DNA Methylation Score     | Calculation of DNA methylation scores            |
| Multiqc                   | Generation of multi-sample quality control report|
+---------------------------+--------------------------------------------------+
""")
workflow.onComplete {
    log.info("""\
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            THE END
                       Have a Great Day!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
""")
}
