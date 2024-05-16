log.info"""\
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  Nextflow Picard DNA Methylation Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *** Steps          
   
   1)  Fastqc
   2)  Trimming fastqc files
   3)  BWA Aligner
   2)  Filter to BAM files                      
   3)  Removal of PCR duplicates                
   4)  Sorting of BAM files                     
   5)  Marking PCR duplicates                   
   6)  Picard tools for processing BAM files    
   7)  Collecting hybrid selection metrics      
   8)  Collecting mark duplicates metrics       
   9)  DNA Methylation Analysis                
   10) Statistics on alignment files
   11) DNA Methylation Matrix
   12) DNAm Full Matrix
   12) Estimate Cell Count 
   13) Picard Scores 
   14) Multiqc          
     
 *** Rules
     
 #Reference and intervals should be ready before using the pipeline.     
   samtools faidx reference.fa
   picard CreateSequenceDictionary REFERENCE=reference.fa OUTPUT=reference.fa.dict
   picard BedToIntervalList 
   I=bed file downloaded from the Twist panel 
   O=covered_targets_Twist_Methylome_reference.intervals 
   SD=reference.dict
       
 *** Parameters Info   
     
     genome_folder       : hg19 (or hg38) BWA Reference folder must be use.
     intervals           : Covered_targets_Twist_Methylome_reference.intervals.
"""

include { Fastqc } from '../modules/Fastqc'
include { Trim_galore } from '../modules/Trim_galore'
include { Alignment } from '../modules/Alignment'
include { Sambamba } from '../modules/Sambamba' 
include { Sorted_Bam_Files } from '../modules/Sorted_Bam_Files' 
include { Mark_duplicated } from '../modules/Mark_duplicated' 
include { Collect_HS_Metrics } from '../modules/Collect_HS_Metrics' 
include { Collect_MM_Metrics } from '../modules/Collect_MM_Metrics' 
include { MethylDackel } from '../modules/MethylDackel'
include { bedGraph } from '../modules/bedGraph'
include { MethylKit } from '../modules/MethylKit'
include { Processed_bedGraph } from '../modules/Processed_bedGraph'
include { Samtools_stats } from '../modules/Samtools_stats'
include { DNAm_Full_Matrix } from '../modules/DNAm_Full_Matrix'
include { Multiqc } from '../modules/Multiqc'
include { Methylation_Matrix } from '../modules/Methylation_Matrix'
include { Estimate_cell_counts } from '../modules/Estimate_cell_counts'
include { DNA_Methylation_Scores } from '../modules/DNA_Methylation_Scores'


workflow Picard_pipeline {

    take:
    reads   
    outdir 
    
    main:
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    Fastqc(read_pairs_ch)
    Trim_galore(read_pairs_ch)
    Alignment(Trim_galore.out.fq)
       myBamSample = Alignment.out.bam
    Sambamba(myBamSample)
       sortedBam = Sambamba.out
    Sorted_Bam_Files(sortedBam)
        sorted_ch = Sorted_Bam_Files.out.sorted_bam
    Mark_duplicated(sorted_ch)
        sorted_mark = Mark_duplicated.out.markdup     
    Collect_HS_Metrics(sorted_mark)
    Collect_MM_Metrics(sorted_mark)
    MethylDackel(sorted_mark)
    bedGraph(sorted_mark)
    bedGraph2 = bedGraph.out
    Processed_bedGraph(bedGraph2)
    Samtools_stats(myBamSample,sorted_mark)
    MethylKit(Mark_duplicated.out.markdup)
    MethylKit2 = MethylKit.out    
    files_ch = MethylKit2.collectFile(name:"*.methylKit", newLine: true)
    Methylation_Matrix(files_ch) 
    DNAm_Full_Matrix(files_ch)
        full_matrix=DNAm_Full_Matrix.out
            full_matrix2=full_matrix.collectFile(name:"*.csv", newLine: true)  
    Estimate_cell_counts(full_matrix2)  
    Meth_Matrix = Methylation_Matrix.out
        files_ch2 = Meth_Matrix.collectFile(name:"*.csv", newLine: true)
    DNA_Methylation_Scores(files_ch2) 
    
    Channel.empty()
          .mix( Fastqc.out )             
          .mix( Trim_galore.out )        
          .mix( Mark_duplicated.out )    
          .mix( Collect_HS_Metrics.out ) 
          .mix( Collect_MM_Metrics.out )     
          .mix( Samtools_stats.out )                                          
          .map { sample_id, files -> files }
          .collect()
          .set { log_files }  

     Multiqc(log_files)
      
}

log.info("""\
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
""")

workflow.onComplete {
    log.info("""\
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            THE END
                       Have a Great Day!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
""")
}
