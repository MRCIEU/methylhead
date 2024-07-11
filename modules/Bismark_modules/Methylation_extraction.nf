process Methylation_extraction {
    
    input:
    tuple val(sample_id), path(bam)
    
    publishDir "${params.outdir}/Methylation" , mode: 'copy'

    output:
    tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz")                 , emit: bedgraph
    tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz")              , emit: coverage
    tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.deduplicated_splitting_report.txt")        , emit: splitting_report
    tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.deduplicated.M-bias.txt")                  , emit: mbias
   
    script:
    """
 bismark_methylation_extractor --comprehensive --bedGraph --gzip --paired-end "${bam}" --CX --no_overlap  
    """
}