process collect_hs_metrics {
 
  publishDir "${params.outdir}/collect-hs-metrics-files/", mode: 'copy'

  input:
    tuple val(sample_id), path(sorted_mark)
    path params.genome_folder
    path interval_file
  
  output:
    tuple val(sample_id), path("${sample_id}_coverage_metrics.txt") , emit: coverage_metrics
    tuple val(sample_id), path("${sample_id}_coverage")             , emit: coverage
         
  script:
    """    
    picard CollectHsMetrics \
    I=${sorted_mark} \
    O=${sample_id}_coverage_metrics.txt \
    R=${params.genome_folder} \
    BAIT_INTERVALS=${interval_file} \
    TARGET_INTERVALS=${interval_file} \
    MINIMUM_MAPPING_QUALITY=20 \
    COVERAGE_CAP=1000 \
    PER_TARGET_COVERAGE=${sample_id}_coverage \
    NEAR_DISTANCE=500
    """     
}


