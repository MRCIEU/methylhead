process Deduplication {
    
    input:
    tuple val(sample_id), path(bam)
    
    publishDir "${params.outdir}/Deduplication" , mode: 'copy'

    output:   
    tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.deduplicated.bam")        , emit: dedup_bam
    tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.deduplication_report.txt"), emit: dedup_report
     
    script:
    """
    deduplicate_bismark --p --bam $bam 
    """
}