
process Picard_mark_duplicated {
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
