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
        files_ch = coverage.collectFile(name:"*.cov.gz", newLine: true)
        script=params.scripts
     Methylation_Matrix(files_ch, script)
     DNAm_Full_Matrix(files_ch, script)
         full_matrix=DNAm_Full_Matrix.out
            full_matrix2=full_matrix.collectFile(name:"*.csv", newLine: true)
     Estimate_cell_counts(full_matrix2, script)    
     Meth_Matrix = Methylation_Matrix.out.meth_matrix
        files_ch2 = Meth_Matrix.collectFile(name:"*.csv", newLine: true)
     DNA_Methylation_Scores(files_ch2, script)
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
          .set { multiqc_files }  
     Multiqc(multiqc_files)

}
