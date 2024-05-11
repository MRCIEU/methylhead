#!/usr/bin/env nextflow

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
     
  Reference and intervals should be ready before using the pipeline.     
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

params.reads = ""
params.genome_folder = ""
params.outdir = "Picard_Results"
params.intervals= "" 

process FASTQC {
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip") , emit: zip
    
    publishDir "${params.outdir}/fastqc/" , mode: 'copy'
    
    script:
    def fastqc_cmd = "fastqc ${reads[0]} ${reads[1]}"  
    "${fastqc_cmd}"

}

process TRIMMING {

    input:
    tuple val(sample_id), path(reads)
  
    publishDir "${params.outdir}/trimmed/" , mode: 'copy'
    
    output:                               
    tuple val(sample_id), path("*_{1,2}.fq.gz") , emit: fq                                  
    tuple val(sample_id), path("*report.txt")                        , emit: log     , optional: true 
 
    script:
    def trim_galore_cmd = "trim_galore --paired ${reads[0]} ${reads[1]} --gzip"
    "${trim_galore_cmd}"

}    

process BWA{

  input:
  tuple val(sample_id), path(fq)
  
  publishDir "${params.outdir}/Bam_files/" , mode: 'copy'
  
  output:
  path("${sample_id}.bam"), emit: bam
      
  script:
  """
  bwameth.py --reference ${params.genome_folder}  ${fq[0]} ${fq[1]} -t 12 | samtools view -b - > ${sample_id}.bam
  """
}

process SAMBAMBA {
    
    input:
    path(myBamSample)
    
    publishDir "${params.outdir}/sambamba/" , mode: 'copy'
    
    output:   
    path("${myBamSample.baseName}_sambamba.bam") , emit: sambamba
    script:
    """
     sambamba view -h -t 16 -T ${params.genome_folder} --filter 'not secondary_alignment and not failed_quality_control and not supplementary and proper_pair and mapping_quality > 0' -f bam -l 0 ${myBamSample} -o ${myBamSample.baseName}_sambamba.bam 
    """
}

process SAMBAMBA_SORT {
    input:
    path bamFile
    
    publishDir "${params.outdir}/sambamba_sorted/" , mode: 'copy'
    
    output:   
    path "${bamFile.baseName}_sorted.bam", emit: sorted
    
    
    script:
    """
    sambamba sort -t 16 -m 30GiB -l 0 ${bamFile} -o ${bamFile.baseName}_sorted.bam 
    """
}

process SAMTOOLS{

    input:
    path sorted_ch
  
    publishDir "${params.outdir}/sambamba/" , mode: 'copy'
    
    output:
    path "${sorted_ch}.bai", emit: bai
    
    script:
    """
    samtools index -@ 16 ${sorted_ch} > ${sorted_ch}.bai
    """
}

process MARK_DUPLICATES {
    input:
    path sortedBam
    
    publishDir "${params.outdir}/picard/", mode: 'copy'
    
    output:   
    path "${sortedBam.baseName}.markdup.bam" , emit: markdup 
    path "${sortedBam.baseName}.picard_markdup_metrics.txt" , emit: picard_markdup_metrics
    
    script:
    """
    picard MarkDuplicates \
    INPUT=${sortedBam} \
    OUTPUT=${sortedBam.baseName}.markdup.bam \
    METRICS_FILE=${sortedBam.baseName}.picard_markdup_metrics.txt \
    REMOVE_DUPLICATES=false \
    ASSUME_SORT_ORDER=coordinate \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500
    """
}


process PICARD_SAMTOOLS {

    input:
    path sorted_mark

    publishDir "${params.outdir}/picard/", mode: 'copy'

    output:
    path "${sorted_mark}.bai", emit: sorted_bai

    script:
    """
    samtools index -@ 16 ${sorted_mark}
    """
}

process PICARD_COLLECT_HSMETRICS {

    input:
    path sorted_mark
    
    output:
    path "${sorted_mark.baseName}_coverage_metrics" , emit: coverage_metrics
    path "${sorted_mark.baseName}_coverage"         , emit: covarage
    
    publishDir "${params.outdir}/picard/CollectHsMetrics/", mode: 'copy'
    
    script:
    """
    picard -Xmx4g -Xms4g CollectHsMetrics \
    I= ${sorted_mark} \
    O= ${sorted_mark.baseName}_coverage_metrics \
    R= ${params.genome_folder} \
    BAIT_INTERVALS= ${params.intervals} \
    TARGET_INTERVALS= ${params.intervals} \
    MINIMUM_MAPPING_QUALITY= 20 \
    COVERAGE_CAP= 1000 \
    PER_TARGET_COVERAGE= ${sorted_mark.baseName}_coverage \
    NEAR_DISTANCE= 500
    """     
}

process PICARD_COLLECT_MMETRICS {
    input:
    path sorted_mark
    
    publishDir "${params.outdir}/picard/CollectMultipleMetrics/", mode: 'copy'
    
    output:
    path "${sorted_mark.baseName}_collectmultiplemetrics.alignment_summary_metrics", emit: alignment_summary_metrics
    path "${sorted_mark.baseName}_collectmultiplemetrics.gc_bias.detail_metrics"   , emit: gc_bias_detail_metrics 
    path "${sorted_mark.baseName}_collectmultiplemetrics.gc_bias.pdf"              , emit: gc_bias
    path "${sorted_mark.baseName}_collectmultiplemetrics.gc_bias.summary_metrics"  , emit: bias_summary_metrics
    path "${sorted_mark.baseName}_collectmultiplemetrics.insert_size_histogram.pdf", emit: size_histogram
    path "${sorted_mark.baseName}_collectmultiplemetrics.insert_size_metrics"      , emit: insert_size_metrics
   
  script:    
   """
   picard -Xmx4g -Xms4g CollectMultipleMetrics \
   I=${sorted_mark} \
   O=${sorted_mark.baseName}_collectmultiplemetrics \
   R=${params.genome_folder} \
   PROGRAM=null \
   PROGRAM=CollectGcBiasMetrics \
   PROGRAM=CollectInsertSizeMetrics \
   PROGRAM=CollectAlignmentSummaryMetrics
    """
}

process MethylDackel {
    input:
    path sorted_mark
    
    publishDir "${params.outdir}/MethylDackel/", mode: 'copy'
    
    output:
    path "${sorted_mark.baseName}.markdup_OT.svg" , emit: markdup_OT
    path "${sorted_mark.baseName}.markdup_OB.svg" , emit: markdup_OB
   

   script:    
   """
   MethylDackel mbias ${params.genome_folder} ${sorted_mark} ${sorted_mark.baseName}.markdup --nOT 0,0,0,98 --nOB 0,0,3,0
    """
}

process MethylDackel_bedGraph{

    input:    
    path sorted_mark
   
    publishDir "${params.outdir}/bedGraph/", mode: 'copy'
    
    output:
    path "${sorted_mark.baseName}_CpG.bedGraph" , emit: bedGrap
    
     script:      
     """
    MethylDackel extract --minDepth 10 --maxVariantFrac 0.25 --nOT 0,0,0,98 --nOB 0,0,3,0 --mergeContext ${params.genome_folder} ${sorted_mark} --keepDupes 
     """

}

process MethylDackel_methylKit{

    input:  
    path sorted_mark
    
    publishDir "${params.outdir}/methylKit/", mode: 'copy'
    
    output:
    path "${sorted_mark.baseName}_CpG.methylKit" , emit: methylKit
    
    script:      
     """
    MethylDackel extract --methylKit ${params.genome_folder} ${sorted_mark} 
     """
}

process Processed_bedGraph {
      
    input:
    path bedGraph

    publishDir "${params.outdir}/Processed_bedGraph/", mode: 'copy'

    output:
    path "${bedGraph}_processed.bedgraph", emit: Processed_bedGraph

    script:
    """
    awk 'BEGIN {FS=OFS=\"\t\"} NR == 1 {print \$0} NR > 1 {print \$1,\$2,\$3,((\$5/(\$5+\$6)*100)+0),\$5,\$6;}' OFMT="%.2f" ${bedGraph} > "${bedGraph}_processed.bedgraph"
    """
}

process Samtools_stats {
    input:
    path(myBamSample)
    path(sorted_mark)
    
    publishDir "${params.outdir}/Samtools_stats/", mode: 'copy'
    
    output:
    path "${myBamSample.baseName}_samtools_stats.txt" , emit: mybam_samtools_stats
    path "${sorted_mark.baseName}_samtools_stats.txt" , emit: sorted_samtools_stats
   
   script:    
   """
   samtools stats ${myBamSample} | grep ^SN | cut -f 2- > ${myBamSample.baseName}_samtools_stats.txt   
   samtools stats ${sorted_mark} | grep ^SN | cut -f 2- > ${sorted_mark.baseName}_samtools_stats.txt
    """
}

process Methylation_Matrix {

    input:
    path methylKit

    output:
    val "picard_methylation.csv", emit: meth_matrix


    script:
    """
    mkdir -p ${baseDir}/${params.outdir}/Methylation_Matrix
    cd ${baseDir}/${params.outdir}/methylKit
    Rscript ${baseDir}/Picard_Methylation_Matrix.R ${baseDir}/${params.outdir}/Methylation_Matrix/
    """
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
  +-----------------------------------------------------------------------------+                
""")
 
workflow {
 read_pairs_ch= Channel.fromFilePairs(params.reads, checkIfExists: true)    
 FASTQC(read_pairs_ch)
 TRIMMING(read_pairs_ch)
 trim_ch=TRIMMING.out.fq
 BWA(trim_ch) 
 myBamSample = BWA.out.bam
 SAMBAMBA(myBamSample)
 sortedBam = SAMBAMBA.out
 SAMBAMBA_SORT(sortedBam) 
 sorted_ch = SAMBAMBA_SORT.out
 SAMTOOLS(sorted_ch) 
 MARK_DUPLICATES(sorted_ch)
 sorted_mark = MARK_DUPLICATES.out.markdup
 PICARD_SAMTOOLS(sorted_mark)
 PICARD_COLLECT_HSMETRICS(sorted_mark)
 PICARD_COLLECT_MMETRICS(sorted_mark)
 MethylDackel(sorted_mark)
 MethylDackel_bedGraph(sorted_mark)
 bedGraph=MethylDackel_bedGraph.out
 Processed_bedGraph(bedGraph)
 MethylDackel_methylKit(sorted_mark)
 methylKit=MethylDackel_methylKit.out 
 files_ch = MethylDackel_methylKit.out.methylKit.collectFile(name:"*.methylKit", newLine: true)
 Methylation_Matrix(files_ch) 
}

workflow.onComplete {
    log.info("""\
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            THE END
                       Have a Great Day!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
""")
}
