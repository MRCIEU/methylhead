process Concordance {
    
    input:
    tuple val(sample_id), path(bam)
    
    publishDir "${params.outdir}/Concordance" , mode: 'copy'

    output:   
    
   tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.deduplicated_all_meth.bam")           ,emit: all_meth_bam
   tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.deduplicated_all_unmeth.bam")         ,emit: all_unmeth_bam
   tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.deduplicated_consistency_report.txt") ,emit: consistency_stats
   tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.deduplicated_mixed_meth.bam")         ,emit: mixed_meth_bam
  
    script:
    """
   methylation_consistency $bam
    """
}



