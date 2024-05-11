
process Picard_collect_hs_metrics {

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