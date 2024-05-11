
process Picard_collect_mm_metrics {
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