#!/usr/bin/env nextflow

include { fastqc } from '../modules/fastqc'
include { trim_galore } from '../modules/trim_galore'
include { interval_file } from '../modules/interval_file'
include { alignment } from '../modules/alignment'
include { sambamba } from '../modules/sambamba' 
include { sorted_bam_files } from '../modules/sorted_bam_files' 
include { mark_duplicated } from '../modules/mark_duplicated' 
include { collect_hs_metrics } from '../modules/collect_hs_metrics' 
include { collect_mm_metrics } from '../modules/collect_mm_metrics' 
include { methyldackel } from '../modules/methyldackel'
include { bedgraph } from '../modules/bedgraph'
include { processed_bedgraph } from '../modules/processed_bedgraph'
include { samtools_stats } from '../modules/samtools_stats'
include { methylkit } from '../modules/methylkit'
include { bsmap_aligment } from '../modules/bsmap_aligment'
include { camda } from '../modules/camda'
include { methylation_matrix_process } from '../modules/methylation_matrix_process'
include { illumina_matrix } from '../modules/illumina_matrix'
include { estimate_cell_counts } from '../modules/estimate_cell_counts'
include { dna_methylation_scores } from '../modules/dna_methylation_scores'
include { camda_matrix } from '../modules/camda_matrix'
include { multiqc } from '../modules/multiqc'
include { qc_report } from '../modules/qc_report'
include { association_test } from '../modules/association_test' 
 
workflow pipeline {

   take:  
   reads   
   outdir 
       
   main:
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    clean_files_ch = read_pairs_ch.filter { sample_id, reads -> 
        def fileSize = reads.collect { it.size() / (1024 * 1024) } 
        fileSize[0] >= 0.500 && fileSize[1] >= 0.500 }
    fastqc(clean_files_ch) 
    trim_galore(clean_files_ch) 
    interval_file(params.panel,params.genome_fasta) 
       genome_fasta=params.genome_fasta
    alignment(trim_galore.out.fq, genome_fasta)
       mybamsample = alignment.out.bam
    sambamba(mybamsample)
       sortedbam = sambamba.out
    sorted_bam_files(sortedbam)
       sorted_ch     = sorted_bam_files.out.sorted_bam
       sorted_ch_bai = sorted_bam_files.out.bai
    mark_duplicated(sorted_ch)
       sorted_mark      = mark_duplicated.out.markdup 
       reference        = params.genome_fasta   
       intervalfile     = interval_file.out.interval_file
       panel            = interval_file.out.panel_bed 
    collect_hs_metrics(sorted_mark, params.genome_fasta, intervalfile)
    collect_mm_metrics(sorted_mark, reference)
    methyldackel(sorted_mark ,reference)
    bedgraph(sorted_mark)
       bedgraph2 = bedgraph.out
    processed_bedgraph(bedgraph2)
    samtools_stats(mybamsample, sorted_mark, sorted_ch, sorted_ch_bai, panel)   
    methylkit(mark_duplicated.out.markdup, genome_fasta)
    bsmap_aligment(trim_galore.out.fq, genome_fasta) 
       bam_camda = bsmap_aligment.out.bam  
    camda(bam_camda)              
       files_ch= methylkit.out.methylKit_CpG
          .map { file -> file.toString() }
          .collectFile(name:"files.csv",newLine:true)
           assembly = params.assembly 
    methylation_matrix_process(files_ch, assembly) 
       methylation_matrix=methylation_matrix_process.out.meth_matrix
    illumina_matrix(methylation_matrix, assembly)   
    estimate_cell_counts(methylation_matrix,params.cell_reference) 
    dna_methylation_scores(methylation_matrix, assembly)                  
       files_ch= camda.out.camda
          .map { file -> file.toString() }
          .collectFile(name:"files.csv",newLine:true)
    camda_matrix(files_ch, assembly)               
          Channel.empty()
          .mix( fastqc.out )             
          .mix( trim_galore.out )        
          .mix( mark_duplicated.out )    
          .mix( collect_hs_metrics.out ) 
          .mix( collect_mm_metrics.out )     
          .mix( samtools_stats.out )                                          
          .map { sample_id, files -> files }
          .collect()
          .set { multiqc_files }  
    multiqc(multiqc_files)   
          reads_ch = multiqc.out.reads
          reads_hs = multiqc.out.reads_hs
          estimate_cell_counts_ch       = estimate_cell_counts.out.estimate_cell_counts
          coverage_matrix_ch            = methylation_matrix_process.out.coverage_matrix
          methylation_matrix_ch         = methylation_matrix_process.out.meth_matrix
          illumina_matrix_ch            = illumina_matrix.out.illumina_matrix 
          camda_ch                      = camda_matrix.out.camda_matrix
          dna_methylation_scores_ch     = dna_methylation_scores.out.dna_methylation_scores
          samtools_read_counts_files_ch = samtools_stats.out.read_counts
            .map { id, f -> f } 
            qc_files_ch = reads_ch  
            .concat(estimate_cell_counts_ch)
            .concat(coverage_matrix_ch)
            .concat(methylation_matrix_ch)
            .concat(illumina_matrix_ch)
            .concat(camda_ch)
            .concat(dna_methylation_scores_ch)
            .concat(reads_hs)
            .concat(samtools_read_counts_files_ch) 
            .map { file -> file.toAbsolutePath().toString() }
            .collectFile(name: "qc_files.csv", newLine: true)
    qc_report(qc_files_ch)
    association_test(qc_files_ch, params.phenotype, params.models)
}



