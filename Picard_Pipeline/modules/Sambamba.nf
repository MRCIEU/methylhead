
process Sambamba {
    
    input:
    tuple val(sample_id),  path(myBamSample)
    
    publishDir "${params.outdir}/Sambamba/" , mode: 'copy'
    
    output:   
    tuple val(sample_id), path("${sample_id}_sambamba.bam") , emit: sambamba
   
    script:
    """
     sambamba view -h -t 16 -T ${params.genome_folder} --filter 'not secondary_alignment and not failed_quality_control and not supplementary and proper_pair and mapping_quality > 0' -f bam -l 0 ${myBamSample} -o ${sample_id}_sambamba.bam
    """

}