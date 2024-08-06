process Collect_HS_Metrics {

    input:
    tuple val(sample_id), path(sorted_mark)
    path params.genome_folder
    
    output:
    tuple val(sample_id), path("${sample_id}_coverage_metrics.txt"), emit: coverage_metrics
    tuple val(sample_id), path("${sample_id}_coverage"), emit: coverage
       
    publishDir "${params.outdir}/CollectHsMetrics/", mode: 'copy'
    
    script:
    """    
    picard CollectHsMetrics \
    I=${sorted_mark} \
    O=${sample_id}_coverage_metrics.txt \
    R=${params.genome_folder} \
    BAIT_INTERVALS=${baseDir}/${params.outdir}/Interval_file/interval_file \
    TARGET_INTERVALS=${baseDir}/${params.outdir}/Interval_file/interval_file \
    MINIMUM_MAPPING_QUALITY=20 \
    COVERAGE_CAP=1000 \
    PER_TARGET_COVERAGE=${sample_id}_coverage \
    NEAR_DISTANCE=500
    """     
}
