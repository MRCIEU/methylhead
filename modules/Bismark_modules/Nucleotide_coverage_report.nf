process Nucleotide_coverage_report {
    
    input:
    tuple val(sample_id), path(bam)
    
    publishDir "${params.outdir}/Nucleotide_coverage_report" , mode: 'copy'

    output:   
    tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.deduplicated.nucleotide_stats.txt")  , emit: nucleotide_stats2
  
    script:
    """
    bam2nuc $bam --genome_folder ${params.genome_folder} 
    """
}




