process collect_hs_metrics {
 
  input:
    tuple val(sample_id), path(sorted_mark)
    path genome_fasta
    path intervalfile
  
  output:
    tuple val(sample_id), path("${sample_id}_coverage_metrics.txt") , emit: coverage_metrics
    tuple val(sample_id), path("${sample_id}_coverage")             , emit: coverage
         
  script:
    """
    FASTA=\$(realpath ${genome_fasta})
 
    picard CollectHsMetrics \
    	   I=${sorted_mark} \
    	   O=${sample_id}_coverage_metrics.txt \
    	   R=\${FASTA} \
    	   BAIT_INTERVALS=${intervalfile} \
    	   TARGET_INTERVALS=${intervalfile} \
    	   MINIMUM_MAPPING_QUALITY=20 \
    	   COVERAGE_CAP=200 \
    	   PER_TARGET_COVERAGE=${sample_id}_coverage \
    	   NEAR_DISTANCE=250
    """     
}


