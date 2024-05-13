process Bismark_alignment {

    input:
    tuple val(sample_id), path(fq)
    val(u_param)
    val(multicore)
    
    output:
    tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_PE_report.txt"), emit: report
    tuple val(sample_id), path("${sample_id}*_val_1_bismark_bt2_pe.nucleotide_stats.txt"), emit: nucstats
   
    publishDir "${params.outdir}/Bismark_Alignment/" , mode: 'copy'
     
    script:
    def bismark_cmd = "bismark --genome_folder ${params.genome_folder} --nucleotide_coverage -1 ${fq[0]} -2 ${fq[1]} --bam" 
    if (u_param) {
        bismark_cmd += " --u ${u_param}"
    }
    if(multicore) {
        bismark_cmd += " --parallel ${multicore}"
    }
    "${bismark_cmd}"
}