
process Sambamba_sort {
    input:
    path bamFile
    
    publishDir "${params.outdir}/sambamba_sorted/" , mode: 'copy'
    
    output:   
    path "${bamFile.baseName}_sorted.bam", emit: sorted
    
    
    script:
    """
    sambamba sort -t 16 -m 30GiB -l 0 ${bamFile} -o ${bamFile.baseName}_sorted.bam 
    """
}
