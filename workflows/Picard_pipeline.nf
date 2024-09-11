#!/usr/bin/env nextflow

include { Fastqc } from '../modules/Picard_modules/Fastqc'
include { Trim_galore } from '../modules/Picard_modules/Trim_galore'
include { Alignment } from '../modules/Picard_modules/Alignment'
include { Sambamba } from '../modules/Picard_modules/Sambamba' 
include { Sorted_Bam_Files } from '../modules/Picard_modules/Sorted_Bam_Files' 
include { Mark_duplicated } from '../modules/Picard_modules/Mark_duplicated' 
include { Interval_file } from '../modules/Picard_modules/Interval_file'
include { Collect_HS_Metrics } from '../modules/Picard_modules/Collect_HS_Metrics' 
include { Collect_MM_Metrics } from '../modules/Picard_modules/Collect_MM_Metrics' 
include { MethylDackel } from '../modules/Picard_modules/MethylDackel'
include { bedGraph } from '../modules/Picard_modules/bedGraph'
include { MethylKit } from '../modules/Picard_modules/MethylKit'
include { Processed_bedGraph } from '../modules/Picard_modules/Processed_bedGraph'
include { Samtools_stats } from '../modules/Picard_modules/Samtools_stats'
include { DNAm_Full_Matrix } from '../modules/Picard_modules/DNAm_Full_Matrix'
include { Multiqc } from '../modules/Picard_modules/Multiqc'
include { Methylation_Matrix } from '../modules/Picard_modules/Methylation_Matrix'
include { Estimate_cell_counts } from '../modules/Picard_modules/Estimate_cell_counts'
include { DNA_Methylation_Scores } from '../modules/Picard_modules/DNA_Methylation_Scores'

workflow Picard_pipeline {

   take:  
   reads   
   outdir 
       
   main:
    
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    Fastqc(read_pairs_ch)
    Trim_galore(read_pairs_ch)   
       genome_folder=params.genome_folder
    Alignment(Trim_galore.out.fq, genome_folder)
       myBamSample = Alignment.out.bam
    Sambamba(myBamSample)
       sortedBam = Sambamba.out
    Sorted_Bam_Files(sortedBam)
       sorted_ch = Sorted_Bam_Files.out.sorted_bam
    Mark_duplicated(sorted_ch)
       sorted_mark = Mark_duplicated.out.markdup 
       reference   = params.genome_folder 
       Interval_file(params.panel,params.genome_folder)
       params.interval_file = Interval_file.out
    Collect_HS_Metrics(sorted_mark, params.genome_folder,params.interval_file)
    Collect_MM_Metrics(sorted_mark, reference)
    MethylDackel(sorted_mark ,reference)
    bedGraph(sorted_mark)
    bedGraph2 = bedGraph.out
    Processed_bedGraph(bedGraph2)
    Samtools_stats(myBamSample,sorted_mark)   
   MethylKit(Mark_duplicated.out.markdup, genome_folder)      
       Channel.empty()
        .mix(MethylKit.out.methylKit_CpG) 
        .map { sample_id, files -> files }                      
        .collect()                                  
        .set { files }
    DNAm_Full_Matrix(files)
       full_matrix=DNAm_Full_Matrix.out
    Methylation_Matrix(full_matrix)   
    Estimate_cell_counts(full_matrix) 
    DNA_Methylation_Scores(full_matrix)      
    Channel.empty()
          .mix( Fastqc.out )             
          .mix( Trim_galore.out )        
          .mix( Mark_duplicated.out )    
          .mix( Collect_HS_Metrics.out ) 
          .mix( Collect_MM_Metrics.out )     
          .mix( Samtools_stats.out )                                          
          .map { sample_id, files -> files }
          .collect()
          .set { multiqc_files }  
    Multiqc(multiqc_files)      
}
