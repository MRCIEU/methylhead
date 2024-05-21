
process Collect_HS_Metrics {

    input:
    tuple val(sample_id) , path (sorted_mark)
    
    output:
    tuple val(sample_id) , path ("${sample_id}_coverage_metrics.txt") , emit: coverage_metrics
    tuple val(sample_id) , path ("${sample_id}_coverage")         , emit: covarage
    
    publishDir "${params.outdir}/CollectHsMetrics/", mode: 'copy'
    
    script:
    """
    picard -Xmx4g -Xms4g CollectHsMetrics \
    I= ${sorted_mark} \
    O= ${sample_id}_coverage_metrics.txt \
    R= ${params.genome_folder} \
    BAIT_INTERVALS= ${params.intervals} \
    TARGET_INTERVALS= ${params.intervals} \
    MINIMUM_MAPPING_QUALITY= 20 \
    COVERAGE_CAP= 1000 \
    PER_TARGET_COVERAGE= ${sample_id}_coverage \
    NEAR_DISTANCE= 500
    """     
}