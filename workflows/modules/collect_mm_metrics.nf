process collect_mm_metrics {

  publishDir "${params.outdir}/collect-multiple-metrics-files/", mode: 'copy'  
  
  input: 
    tuple val(sample_id) , path (sorted_mark)
    path reference    

  output:
    tuple val(sample_id) , path ("${sample_id}_collectmultiplemetrics.alignment_summary_metrics")    , emit: alignment_summary_metrics
    tuple val(sample_id) , path ("${sample_id}_collectmultiplemetrics.gc_bias.detail_metrics")       , emit: gc_bias_detail_metrics 
    tuple val(sample_id) , path ("${sample_id}_collectmultiplemetrics.gc_bias.pdf")                  , emit: gc_bias
    tuple val(sample_id) , path ("${sample_id}_collectmultiplemetrics.gc_bias.summary_metrics")      , emit: bias_summary_metrics
    tuple val(sample_id) , path ("${sample_id}_collectmultiplemetrics.insert_size_histogram.pdf")    , emit: size_histogram
    tuple val(sample_id) , path ("${sample_id}_collectmultiplemetrics.insert_size_metrics")          , emit: insert_size_metrics
   
  script:    
   """
    picard CollectMultipleMetrics \
    I=${sorted_mark} \
    O=${sample_id}_collectmultiplemetrics \
    R=${reference} \
    PROGRAM=null \
    PROGRAM=CollectGcBiasMetrics \
    PROGRAM=CollectInsertSizeMetrics \
    PROGRAM=CollectAlignmentSummaryMetrics
    """
}
