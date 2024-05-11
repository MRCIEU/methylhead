
process Sambamba {
    
    input:
    path(myBamSample)
    
    publishDir "${params.outdir}/sambamba/" , mode: 'copy'
    
    output:   
    path("${myBamSample.baseName}_sambamba.bam") , emit: sambamba
    script:
    """
     sambamba view -h -t 16 -T ${params.genome_folder} --filter 'not secondary_alignment and not failed_quality_control and not supplementary and proper_pair and mapping_quality > 0' -f bam -l 0 ${myBamSample} -o ${myBamSample.baseName}_sambamba.bam 
    """
}