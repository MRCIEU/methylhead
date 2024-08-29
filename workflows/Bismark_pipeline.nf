#!/usr/bin/env nextflow

include { Fastqc } from '../modules/Bismark_modules/Fastqc'
include { Trim_galore } from '../modules/Bismark_modules/Trim_galore'
include { Alignment } from '../modules/Bismark_modules/Alignment'
include { Deduplication } from '../modules/Bismark_modules/Deduplication'
include { Methylation_extraction } from '../modules/Bismark_modules/Methylation_extraction'
include { DNAm_Full_Matrix } from '../modules/Bismark_modules/DNAm_Full_Matrix'
include { Methylation_Matrix } from '../modules/Bismark_modules/Methylation_Matrix'
include { Estimate_cell_counts } from '../modules/Bismark_modules/Estimate_cell_counts'
include { DNA_Methylation_Scores } from '../modules/Bismark_modules/DNA_Methylation_Scores'
include { Reports } from '../modules/Bismark_modules/Reports'
include { Multiqc } from '../modules/Bismark_modules/Multiqc'

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
 	   genome_folder=params.genome_folder
     bam_files_ch = Alignment(trim_ch,u_param,multicore, genome_folder)
        ch_bam = Alignment.out.alignment_bam
     Deduplication(ch_bam) 
         dedup_bam=Deduplication.out.dedup_bam
     Methylation_extraction(dedup_bam)
     coverage= Methylation_extraction.out.coverage               
       Channel.empty()
         .mix( Methylation_extraction.out )
         .map { sample_id, files -> files }
	 .collect()
	 .set{ R_files }     
    Methylation_Matrix(R_files)
    DNAm_Full_Matrix(R_files)
    Estimate_cell_counts(R_files)    
    DNA_Methylation_Scores(fR_files)
     	Channel.empty()
	     .mix( Alignment.out.alignment_report )
             .mix( Deduplication.out.dedup_report )
             .mix( Methylation_extraction.out.splitting_report )
             .mix( Methylation_extraction.out.mbias )
  	     .map { sample_id, files -> files}
             .collect()
             .set { Report_files }
     Reports(Report_files) 
        Channel.empty()
          .mix( Fastqc.out )
          .mix( Trim_galore.out )
          .mix( Alignment.out )
          .mix( Deduplication.out )
          .mix( Methylation_extraction.out ) 
          .map { sample_id, files -> files }
          .collect()
          .set { multiqc_files }  
     Multiqc(multiqc_files)
}
