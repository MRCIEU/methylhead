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
include { Bwa_aligner } from '../modules/Bwa_aligner'
include { Sambamba } from '../modules/Sambamba' 
include { Sambamba_sort } from '../modules/Sambamba_sort' 
include { Samtools } from '../modules/Samtools' 
include { Picard_mark_duplicated } from '../modules/Picard_mark_duplicated' 
include { Picard_samtools } from '../modules/Picard_samtools' 
include { Picard_collect_hs_metrics } from '../modules/Picard_collect_hs_metrics' 
include { Picard_collect_mm_metrics } from '../modules/Picard_collect_mm_metrics' 
include { MethylDackel } from '../modules/MethylDackel'
include { MethylDackel_bedGraph } from '../modules/MethylDackel_bedGraph'
include { MethylDackel_methylKit } from '../modules/MethylDackel_methylKit'
include { Processed_bedGraph } from '../modules/Processed_bedGraph'
include { Samtools_stats } from '../modules/Samtools_stats'
include { Methylation_Matrix } from '../modules/Methylation_Matrix'


workflow Picard_pipeline {

    take:
    reads   
    outdir
    
    main:
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    Fastqc(read_pairs_ch)
    Trim_galore(read_pairs_ch)
    Bwa_aligner(Trim_galore.out.fq)
       myBamSample = Bwa_aligner.out.bam
    Sambamba(myBamSample)
       sortedBam = Sambamba.out
    Sambamba_sort(sortedBam)
        sorted_ch = Sambamba_sort.out
    Samtools(sorted_ch) 
    Picard_mark_duplicated(sorted_ch)
        sorted_mark = Picard_mark_duplicated.out.markdup    
    Picard_samtools(sorted_mark)   
    Picard_collect_hs_metrics(sorted_mark)
    Picard_collect_mm_metrics(sorted_mark)
    MethylDackel(sorted_mark)
    MethylDackel_bedGraph(sorted_mark)
    bedGraph = MethylDackel_bedGraph.out
    Processed_bedGraph(bedGraph)
    MethylDackel_methylKit(Picard_mark_duplicated.out.markdup)
    MethylKit = MethylDackel_methylKit.out 
    files_ch = MethylKit.collectFile(name:"*.methylKit", newLine: true)
    Methylation_Matrix(files_ch)
    Estimate_cell_counts(files_ch)
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
  |Methylation Matrix       |Creating Methylation Matrix with R script          |
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
