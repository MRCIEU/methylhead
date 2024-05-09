#!/usr/bin/env nextflow

log.info"""\
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  Nextflow DNA Methylation Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *** Steps          
     1. QC Steps
     2. Bismark
     3. Methylation Matrix  
     4. DNA methylation indices of exposure and phenotype
 *** Rules    
     All the process names should be written with BIG LETTERS.
     All workflows are processed sequentially; they cannot be used separately.       
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
params.reads = ""
params.genome_folder = ""
params.outdir = "Bismark_Results"
params.u_param = "" 
params.t_param = ""
params.memory_param = "" 
params.multicore= ""
params.cores= ""

process FASTQC {
      
    input:
    tuple val(sample_id), path(reads)
    val(t_param)
    val(memory_param)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip") , emit: zip
    
    publishDir "${params.outdir}/fastqc/" , mode: 'copy'
    
    script:
    def fastqc_cmd = "fastqc ${reads[0]} ${reads[1]} --output_dir ${params.outdir}"
    if (t_param) {
        fastqc_cmd += " --threads ${params.t_param}"
    }
    if (memory_param) {
        fastqc_cmd += " --memory ${params.memory_param}"
    }    
    "${fastqc_cmd}"

}

process TRIMMING {
    input:
    tuple val(sample_id), path(reads)
    val(cores)
    
    output:                               
    tuple val(sample_id), path("*_{1,2}.fq.gz") , emit: fq                                  
    tuple val(sample_id), path("*report.txt")                        , emit: log     , optional: true
    tuple val(sample_id), path("*.html")                             , emit: html    , optional: true   
   
    publishDir "${params.outdir}/trimmed/" , mode: 'copy'
    
 
    script:
    def trim_galore_cmd = "trim_galore --paired ${reads[0]} ${reads[1]} --gzip"
    if(cores) {
       trim_galore_cmd += " --cores ${cores}"
    }
    "${trim_galore_cmd}"

}

process BISMARK {

    input:
    tuple val(sample_id), path(fq)
    val(u_param)
    val(multicore)
    
    output:
    tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_PE_report.txt"), emit: report
    tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.nucleotide_stats.txt"), emit: nucstats
   
    publishDir "${params.outdir}/Bismark_Alignment/" , mode: 'copy'
     
    script:
    def bismark_cmd = "bismark --genome_folder ${params.genome_folder} --nucleotide_coverage -1 ${fq[0]} -2 ${fq[1]} --bam" 
    if (u_param) {
        bismark_cmd += " --u ${u_param}"
    }
    if(multicore) {
        bismark_cmd += " --parallel ${multicore}"
    }
    "${bismark_cmd}"
}
process DEDUPLICATION {
    
    input:
    tuple val(sample_id), path(bam)
    
    publishDir "${params.outdir}/Deduplication" , mode: 'copy'

    output:   
    tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.deduplicated.bam")        , emit: bam
    tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.deduplication_report.txt"), emit: report
     
    script:
    """
    deduplicate_bismark --p --bam $bam 
    """
}
process METHYLATION {
    
    input:
    tuple val(sample_id), path(bam)
    
    publishDir "${params.outdir}/Methylation" , mode: 'copy'

    output:
    tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz")                 , emit: bedgraph
    tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz")              , emit: coverage
    tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.deduplicated_splitting_report.txt")        , emit: report
    tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.deduplicated.M-bias.txt")                  , emit: mbias
   
    script:
    """
 bismark_methylation_extractor --comprehensive --bedGraph --gzip --paired-end "${bam}" --CX --no_overlap  
    """
}

log.info("""\
  +---------------------------------------------------------------+
  | Pipeline Step    | Description                                |
  +------------------+--------------------------------------------+
  | Fastqc           | Fastq Quality Control                      |
  | Trimming         | Trimming low-quality bases                 |
  |---------------------------------------------------------------| 
  | Alignment        | Paired-end fastq alignment to bam files    |
  | Deduplication    | Removal of PCR duplicates                  |
  | DNA Methylation  | DNA Methylation Analysis                   |
  | Report           | Bismark Processing Report                  |
  |---------------------------------------------------------------|
  |                   Methylation Matrix                          |
  |---------------------------------------------------------------|
  |    DNA methylation indices of exposure and phenotype          |
  +---------------------------------------------------------------+
""")

workflow {

 read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)    
 t_param=params.t_param
 memory_param=params.memory_param
 FASTQC(read_pairs_ch,params.t_param,params.memory_param)
 cores=params.cores
 TRIMMING(read_pairs_ch,cores)
 trim_ch=TRIMMING.out.fq
 u_param= params.u_param 
 multicore= params.multicore 
 bam_files_ch = BISMARK(trim_ch,params.u_param,multicore)
 ch_bam = BISMARK.out.bam
 DEDUPLICATION(ch_bam) 
 dedup_bam=DEDUPLICATION.out.bam
 METHYLATION(dedup_bam)
 coverage= METHYLATION.out                      
}

workflow.onComplete {
    log.info("""\
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            THE END
                       Have a Great Day!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
""")
}
