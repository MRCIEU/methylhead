
process Samtools {

    input:
    path sorted_ch
  
    publishDir "${params.outdir}/sambamba/" , mode: 'copy'
    
    output:
    path "${sorted_ch}.bai", emit: bai
    
    script:
    """
    samtools index -@ 16 ${sorted_ch} > ${sorted_ch}.bai
    """
}
