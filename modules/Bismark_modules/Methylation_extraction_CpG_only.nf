process Methylation_extraction_CpG_only {
    
    input:
    tuple val(sample_id), path(bam)
    
    publishDir "${params.outdir}/Methylation_CpG_only" , mode: 'copy'

    output:
    
    tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz")              , emit: coverage
  
    script:
    """
    bismark_methylation_extractor --comprehensive --bedGraph --gzip --paired-end "${bam}" --no_overlap  
    
    """
}