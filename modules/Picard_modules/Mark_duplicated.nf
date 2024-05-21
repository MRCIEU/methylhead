
process Mark_duplicated {
    
    input:
    
    tuple val(sample_id),  path (sortedBam)
    
    publishDir "${params.outdir}/Mark_duplicated/", mode: 'copy'
    
    output:   
    tuple val(sample_id) , path("${sample_id}.markdup.bam") , emit: markdup 
    tuple val(sample_id) , path("${sample_id}.picard_markdup_metrics.txt") , emit: picard_markdup_metrics
    
    script:
    """
    picard MarkDuplicates \
    INPUT=${sortedBam} \
    OUTPUT=${sample_id}.markdup.bam \
    METRICS_FILE=${sample_id}.picard_markdup_metrics.txt \
    REMOVE_DUPLICATES=false \
    ASSUME_SORT_ORDER=coordinate \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500
            
    samtools index -@ 16 ${sample_id}.markdup.bam  > ${sample_id}.markdup.bam.bai
    """
}
