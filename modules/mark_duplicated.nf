process mark_duplicated {
  
  input: 
    tuple val(sample_id),  path (sortedbam)
    
  output:   
    tuple val(sample_id) , path("${sample_id}.markdup.bam")                , emit: bam 
    tuple val(sample_id) , path("${sample_id}.markdup_metrics.txt")        , emit: metrics
    tuple val(sample_id) , path("${sample_id}.markdup.bam.bai")            , emit: bai
    
  script:
    """
    picard MarkDuplicates \
    INPUT=${sortedbam} \
    OUTPUT=${sample_id}.markdup.bam \
    METRICS_FILE=${sample_id}.markdup_metrics.txt \
    REMOVE_DUPLICATES=false \
    ASSUME_SORT_ORDER=coordinate \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500
          
    samtools index -@ 16 ${sample_id}.markdup.bam > ${sample_id}.markdup.bam.bai
    """
}
