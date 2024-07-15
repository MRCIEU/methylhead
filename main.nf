#!/usr/bin/env nextflow

params.pipeline = ''

if (params.pipeline != 'picard' && params.pipeline != 'bismark') {
    error "OppSs!..wait wait wait! We just have two pipelines: Please choose either 'picard' or 'bismark'"
}

log.info "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
log.info "                  Nextflow DNA Methylation Pipeline"
log.info "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

if (params.pipeline == 'bismark') {
    log.info """
     *** Bismark Pipeline Steps
     1. QC Steps
     2. Alignment & Methylation steps
     3. Methylation Matrix
     4. Concordance   
     5. DNAm Full Matrix 
     6. Estimation cell counts
     7. DNA Methylation Score
     8. Bismark Reports
     9. Multiqc 
    """
} else if (params.pipeline == 'picard') {
    log.info """
     *** Picard Pipeline Steps      
     1.  Fastqc
     2   Trimming fastqc files
     3.  BWA Aligner
     4.  Filter to BAM files                      
     5.  Removal of PCR duplicates                
     6.  Sorting of BAM files                     
     7.  Marking PCR duplicates                   
     8.  Picard tools for processing BAM files    
     9.  Collecting hybrid selection metrics      
     10. Collecting mark duplicates metrics       
     11. DNA Methylation Analysis                
     12. Statistics on alignment files
     13. DNA Methylation Matrix
     14. DNAm Full Matrix
     15. Estimate Cell Count 
     16. DNA Methylation Score
     17. Multiqc 
    """
}

include { Bismark_pipeline } from './workflows/Bismark_pipeline'
include { Picard_pipeline } from './workflows/Picard_pipeline'

workflow {
    reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
    outdir = params.outdir

    if (params.pipeline == 'bismark') {
        Bismark_pipeline(reads, outdir)
    } else if (params.pipeline == 'picard') {
        Picard_pipeline(reads, outdir)
    } else {
        error "Invalid pipeline selected. Choose either 'bismark' or 'picard'."
    }
}

workflow.onComplete {
    if (params.pipeline == 'bismark') {
log.info("""
          +---------------------------+--------------------------------------------------+
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
          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                    Bismark Pipeline Completed                     
                                       Have a Great Day!                          
          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~         
""")     

} else if (params.pipeline == 'picard') {

log.info("""
           +-----------------------------------------------------------------------------+
           |      Pipeline Step      |                  Description                      |
           |-----------------------------------------------------------------------------| 
           |FastQC                   |Quality control for raw sequencing data            |
           |Trimming FastQC files    |Quality control for trimmed sequencing data        |
           |BWA Aligner              |Alignment of sequencing reads to a reference genome|
           |Filter to BAM files      |Conversion of aligned reads to BAM format          |
           |Sambamba                 |Filter to bam files                                |
           |Deduplication            |Removal of PCR duplicates                          |
           |Samtools                 |Sorting of bam files                               |
           |Mark Duplicates          |Marking PCR duplicates                             |
           |Picard Samtools          |Picard tools for processing bam files              |
           |Picard Collect HS Metrics|Collecting hybrid selection metrics                |
           |Picard Collect MM Metrics|Collecting mark duplicates metrics                 |
           |MethylDackel             |DNA Methylation Analysis                           |
           |MethylDackel bedGraph    |Produce bedGraph for MethylDackel                  |
           |MethylDackel methylKit   |Produce cytosine report for methylKit              |
           |Processed_bedGraph       |Processed bedGraph                                 |
           |Samtools Stats           |Statistics on alignment files                      |
           |Methylation Matrix       |Creating Methylation Matrix                        |
           |DNAm Full Matrix         |Creating Full DNAm Matrix                          |
           |Estimate Cell Count      |Estimate Cell Counts from Loyfer Data              |
           |DNA Methylation Scores   |Creating Metyhlation Scores with R script          | 
           |Multiqc                  |Creating Multi QC report                           |
           +-----------------------------------------------------------------------------+              
           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                    Picard Pipeline Completed                    
                                       Have a Great Day!                          
           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~              
""")
 }
}


workflow.onComplete {

    def msg = """\

        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()
}
