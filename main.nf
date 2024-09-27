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
       1. **QC Steps**  
          - Fastqc: Quality control for raw sequencing data  
          - Trim_galore: Trimming low-quality bases from raw sequencing reads  
       2. **Alignment & Methylation steps**  
          - Alignment: Aligning sequencing reads using Bismark  
          - Deduplication: Removing PCR duplicates  
          - Methylation_extraction: Extracting methylation information  
          - BSmap_Aligment: Aligning bisulfite sequencing data  
       3. **Concordance**  
          - Concordance score from Bismark tool
          - CAMDA: Processing CAMDA dataset  
       4. **DNAm Full Matrix**  
          - DNAm_Matrix: Creating full DNA methylation matrix  
          - Illumina_Matrix: Creating Illumina-based matrix  
       5. **Estimation of cell counts**  
          - Estimate_cell_counts: Estimating cell counts using methylation data  
       6. **DNA Methylation Score**  
          - DNA_Methylation_Scores: Generating methylation scores using R  
       7. **Bismark Reports**  
          - Nucleotide coverage report from Bismark tool
          - Reports: Generating reports for alignment and methylation steps  
       8. **MultiQC**  
          - Multiqc: Aggregating and summarizing results into a final MultiQC report  
    """
} else if (params.pipeline == 'picard') {
    log.info """
      *** Picard Pipeline Steps      
       1. **FastQC**  
          - Quality control for raw sequencing data.  
       2. **Trim Galore**  
          - Trimming low-quality bases from raw sequencing reads.  
       3. **Interval file**  
          - Creating an interval file for targeted sequencing regions.  
       4. **bwa_meth Aligner**  
          - Aligning sequencing reads to the reference genome using BWA-Meth.  
       5. **Sambamba for filtering to BAM files**  
          - Filtering aligned reads and converting to BAM format.  
       6. **Picard: Sorting BAM Files**  
          - Sorting BAM files by coordinates using Picard tools.  
       7. **Picard: Removing PCR duplicates**  
          - Removing PCR duplicates from BAM files using Picard.  
       8. **Picard: Collect HS Metrics (QC)**  
          - Collecting hybrid selection metrics for quality control.  
       9. **Picard: Collect MM Metrics (QC)**  
          - Collecting metrics related to mark duplicates for quality control.  
       10. **MethylDackel (QC)**  
           - DNA methylation analysis and quality control.  
       11. **bedGraph**  
           - Producing bedGraph files for visualizing methylation data.  
       12. **Processed bedGraph**  
           - Processing bedGraph files for further analysis.  
       13. **Samtools Stats**  
           - Generating statistics on alignment files using Samtools.  
       14. **MethylKit**  
           - Producing methylation reports from cytosine sites using MethylKit.  
       15. **BSmap Alignment for CAMDA**  
           - Aligning bisulfite sequencing data for CAMDA analysis.  
       16. **CAMDA**  
           - Processing CAMDA dataset for concordance analysis.  
       17. **DNAm Matrix**  
           - Creating the full DNA methylation matrix.  
       18. **Illumina Matrix (450k annotation)**  
           - Creating a matrix based on Illumina 450k annotation data.  
       19. **Estimate Cell Counts**  
           - Estimating cell counts for CD4T, CD8T, NK, Mono, B cells, and Granulocytes.  
       20. **DNA Methylation Scores**  
           - Generating methylation scores using R for further analysis.  
       21. **Camda Matrix**  
           - Creating a matrix for CAMDA analysis results.  
       22. **MultiQC**  
           - Aggregating and summarizing results into a comprehensive MultiQC report.
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
       

           +--------------------------------------------------------------------+
           | Pipeline Step         | Description                                |
           |--------------------------------------------------------------------|
           | FastQC                | Quality control for raw sequencing data    |                            
           | Trim Galore           | Trimming low-quality bases from raw seq    |
           | Alignment             | Alignment of sequencing reads using Bismark|                               
           | Deduplication         | Removal of PCR duplicates                  |
           | Methylation Extraction| Extracting methylation information         |
           | BSmap Alignment       | Aligning bisulfite sequencing data         |
           | CAMDA                 | Processing CAMDA dataset using python      |
           | DNAm Matrix           | Creating full DNA methylation matrix       |
           | Illumina Matrix       | Creating methylation matrix using R        |
           | Estimate Cell Counts  | Estimating cell counts using R             |                   
           | DNA Methylation Scores| Generating methylation scores using R      |
           | Camda Matrix          | Creating matrix for CAMDA analysis using R |
           | Reports               | Generating reports for analysis steps      |
           | MultiQC               | Aggregating results into a MultiQC Report  |
           +--------------------------------------------------------------------+

           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                  Bismark Pipeline Completed
                                       Have a Great Day!
           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     
""")     

} else if (params.pipeline == 'picard') {

log.info("""
           +--------------------------------------------------------------------------+
           | Pipeline Step         | Description                                      |
           |--------------------------------------------------------------------------|
           | FastQC                | Quality control for raw sequencing data          |
           | Trim Galore           | Trimming low-quality bases from raw reads        |
           | Interval file         | Interval file creation for targeted seq          |
           | Alignment             | Alignment of sequencing reads to a ref file      |
           | Sambamba              | Sorting and filtering of BAM files               |
           | Sorted BAM Files      | Sorting of BAM files using Samtools              |
           | Mark Duplicates       | Marking PCR duplicates in BAM files              |
           | Collect HS Metrics    | Collecting hybrid selection metrics              |
           | Collect MM Metrics    | Collecting mark duplicates metrics               |
           | MethylDackel          | DNA methylation analysis                         |
           | bedGraph              | Producing bedGraph files for methylation         |
           | Processed bedGraph    | Further processing of bedGraph files             |
           | Samtools Stats        | Generating statistics for alignment              |
           | MethylKit             | Producing methylation reports from cytosine sites|
           | BSmap Alignment       | Aligning bisulfite sequencing data               |
           | CAMDA                 | Processing CAMDA dataset using python            |
           | DNAm Matrix           | Creating DNA methylation matrix using R          |
           | Illumina Matrix       | Creating matrix from Illumina data using R       |
           | Estimate Cell Counts  | Estimating cell counts using R                   |
           | DNA Methylation Scores| Generating methylation scores using R            |
           | Camda Matrix          | Creating matrix for CAMDA analysis using R       |
           | MultiQC               | Aggregating results into a multiQC report        |
           +--------------------------------------------------------------------------+
           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                    Picard Pipeline Completed
                                       Have a Great Day!
           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
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
