#!/usr/bin/env nextflow

include { Fastqc } from '../modules/Picard_modules/Fastqc'
include { Trim_galore } from '../modules/Picard_modules/Trim_galore'
include { Interval_file } from '../modules/Picard_modules/Interval_file'
include { Alignment } from '../modules/Picard_modules/Alignment'
include { Sambamba } from '../modules/Picard_modules/Sambamba' 
include { Sorted_Bam_Files } from '../modules/Picard_modules/Sorted_Bam_Files' 
include { Mark_duplicated } from '../modules/Picard_modules/Mark_duplicated' 
include { Collect_HS_Metrics } from '../modules/Picard_modules/Collect_HS_Metrics' 
include { Collect_MM_Metrics } from '../modules/Picard_modules/Collect_MM_Metrics' 
include { MethylDackel } from '../modules/Picard_modules/MethylDackel'
include { bedGraph } from '../modules/Picard_modules/bedGraph'
include { Processed_bedGraph } from '../modules/Picard_modules/Processed_bedGraph'
include { Samtools_stats } from '../modules/Picard_modules/Samtools_stats'
include { MethylKit } from '../modules/Picard_modules/MethylKit'
include { BSmap_Aligment } from '../modules/Picard_modules/BSmap_Aligment'
include { CAMDA } from '../modules/Picard_modules/CAMDA'
include { DNAm_Matrix } from '../modules/Picard_modules/DNAm_Matrix'
include { Illumina_Matrix } from '../modules/Picard_modules/Illumina_Matrix'
include { Estimate_cell_counts } from '../modules/Picard_modules/Estimate_cell_counts'
include { DNA_Methylation_Scores } from '../modules/Picard_modules/DNA_Methylation_Scores'
include { Camda_matrix } from '../modules/Picard_modules/Camda_matrix'
include { Multiqc } from '../modules/Picard_modules/Multiqc'
include { QC_Report } from '../modules/Picard_modules/QC_Report'
include { Association_test } from '../modules/Picard_modules/Association_test' 
 

workflow Picard_pipeline {

   take:  
   reads   
   outdir 
       
   main:
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    clean_files_ch = read_pairs_ch.filter { sample_id, reads -> 
        def fileSize = reads.collect { it.size() / (1024 * 1024) } 
        fileSize[0] >= 0.500 && fileSize[1] >= 0.500 }
    Fastqc(clean_files_ch) 
    Trim_galore(clean_files_ch) 
    Interval_file(params.panel,params.genome_folder) 
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
       interval_file = Interval_file.out
    Collect_HS_Metrics(sorted_mark, params.genome_folder,interval_file)
    Collect_MM_Metrics(sorted_mark, reference)
    MethylDackel(sorted_mark ,reference)
    bedGraph(sorted_mark)
    bedGraph2 = bedGraph.out
    Processed_bedGraph(bedGraph2)
    Samtools_stats(myBamSample,sorted_mark)   
    MethylKit(Mark_duplicated.out.markdup, genome_folder)
    BSmap_Aligment(Trim_galore.out.fq, genome_folder) 
       bam_camda = BSmap_Aligment.out.bam  
    CAMDA(bam_camda)              
       files_ch= MethylKit.out.methylKit_CpG
          .map { file -> file.toString() }
          .collectFile(name:"files.csv",newLine:true)
    DNAm_Matrix(files_ch) 
       full_matrix=DNAm_Matrix.out.meth_matrix
    Illumina_Matrix(full_matrix)   
    Estimate_cell_counts(full_matrix) 
    DNA_Methylation_Scores(full_matrix)                  
       files_ch= CAMDA.out.camda
          .map { file -> file.toString() }
          .collectFile(name:"files.csv",newLine:true)
    Camda_matrix(files_ch)               
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
          reads_ch = Multiqc.out.reads
          estimate_cell_counts_ch = Estimate_cell_counts.out.estimate_cell_counts
          coverage_matrix_ch = DNAm_Matrix.out.coverage_matrix
          methylation_matrix_ch = DNAm_Matrix.out.meth_matrix
          Illumina_Matrix_ch = Illumina_Matrix.out.Illumina_matrix  
          camda_ch = Camda_matrix.out.camda_matrix
          dna_methylation_scores_ch = DNA_Methylation_Scores.out.dna_methylation_scores
            qc_files_ch = reads_ch
            .concat(estimate_cell_counts_ch)
            .concat(coverage_matrix_ch)
            .concat(methylation_matrix_ch)
            .concat(Illumina_Matrix_ch)
            .concat(camda_ch)
            .concat(dna_methylation_scores_ch)
            .map { file -> file.toAbsolutePath().toString() }
            .collectFile(name: "qc_files.csv", newLine: true)
    QC_Report(qc_files_ch)
   Association_test(qc_files_ch,params.phenotype, params.models)
}



