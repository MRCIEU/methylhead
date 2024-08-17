
process Sorted_Bam_Files {

    input:
    
    tuple val(sample_id),  path (sortedBam)
    
    publishDir "${params.outdir}/Sorted_Bam_Files/" , mode: 'copy'
    
    output:   
   
    tuple val(sample_id) , path ("${sample_id}_sorted.bam"), emit: sorted_bam
    tuple val(sample_id) , path ("${sample_id}_sorted.bam.bai"), emit: bai
    
    
    script:
    """
    sambamba sort -t 16 -m 16GiB -l 0 ${sortedBam} -o ${sample_id}_sorted.bam    
    samtools index -@ 16 ${sample_id}_sorted.bam  > ${sample_id}_sorted.bam.bai
    """
}
