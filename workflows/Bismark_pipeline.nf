#!/usr/bin/env nextflow

include { Fastqc } from '../modules/Bismark_modules/Fastqc'
include { Trim_galore } from '../modules/Bismark_modules/Trim_galore'
include { Alignment } from '../modules/Bismark_modules/Alignment'
include { Deduplication } from '../modules/Bismark_modules/Deduplication'
include { Methylation_extraction } from '../modules/Bismark_modules/Methylation_extraction'
include { Concordance } from '../modules/Bismark_modules/Concordance'
include { BSmap_Aligment } from '../modules/Bismark_modules/BSmap_Aligment'
include { CAMDA } from '../modules/Bismark_modules/CAMDA'
include { DNAm_Matrix } from '../modules/Bismark_modules/DNAm_Matrix'
include { Illumina_Matrix } from '../modules/Bismark_modules/Illumina_Matrix'
include { Estimate_cell_counts } from '../modules/Bismark_modules/Estimate_cell_counts'
include { DNA_Methylation_Scores } from '../modules/Bismark_modules/DNA_Methylation_Scores'
include { Camda_matrix } from '../modules/Bismark_modules/Camda_matrix'
include { Reports } from '../modules/Bismark_modules/Reports'
include { Nucleotide_coverage_report } from '../modules/Bismark_modules/Nucleotide_coverage_report'
include { Multiqc } from '../modules/Bismark_modules/Multiqc'

workflow Bismark_pipeline {
    take:
    reads   
    outdir 
    
    main:
    
     read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
     clean_files_ch = read_pairs_ch.filter { sample_id, reads -> 
        def fileSize = reads.collect { it.size() / (1024 * 1024) } 
        fileSize[0] >= 100 && fileSize[1] >= 100  }
     Fastqc(clean_files_ch) 
         cores=params.cores
     Trim_galore(clean_files_ch,cores)
         trim_ch=Trim_galore.out.fq         
         u_param= params.u_param 
         multicore= params.multicore
 	   genome_folder=params.genome_folder
     bam_files_ch = Alignment(trim_ch,u_param,multicore, genome_folder)
         ch_bam = Alignment.out.alignment_bam
     Deduplication(ch_bam) 
         dedup_bam=Deduplication.out.dedup_bam
     Methylation_extraction(dedup_bam)
     Concordance(dedup_bam)
     BSmap_Aligment(trim_ch, genome_folder)
         bam_camda = BSmap_Aligment.out.bam
     CAMDA(bam_camda) 
         files_ch= Methylation_extraction.out.coverage
             .map { file -> file.toString() }
             .collectFile(name:"files.csv",newLine:true)
     DNAm_Matrix(files_ch) 
         full_matrix=DNAm_Matrix.out.meth_matrix
     Illumina_Matrix(full_matrix)   
     Estimate_cell_counts(full_matrix) 
     DNA_Methylation_Scores(full_matrix)                        
         files_ch_camda= CAMDA.out.camda
             .map { file -> file.toString() }
             .collectFile(name:"files.csv",newLine:true)
     Camda_matrix(files_ch_camda)   
    	   Channel.empty()
	           .mix( Alignment.out.alignment_report )
             .mix( Deduplication.out.dedup_report )
             .mix( Methylation_extraction.out.splitting_report )
             .mix( Methylation_extraction.out.mbias )
  	         .map { sample_id, files -> files}
             .collect()
             .set { Report_files }
     Nucleotide_coverage_report(dedup_bam)
     Reports(Report_files)             
          Channel.empty()
             .mix( Fastqc.out )
             .mix( Trim_galore.out )
             .mix( Alignment.out )
             .mix( Deduplication.out )
             .mix( Methylation_extraction.out ) 
             .mix( Concordance.out )
             .mix( Nucleotide_coverage_report.out ) 
             .map { sample_id, files -> files }
             .collect()
             .set { multiqc_files }     
     Multiqc(multiqc_files)
}
