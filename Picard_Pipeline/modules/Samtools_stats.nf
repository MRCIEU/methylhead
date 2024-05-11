process Samtools_stats {
    input:
    path(myBamSample)
    path(sorted_mark)
    
    publishDir "${params.outdir}/Samtools_stats/", mode: 'copy'
    
    output:
    path "${myBamSample.baseName}_samtools_stats.txt" , emit: mybam_samtools_stats
    path "${sorted_mark.baseName}_samtools_stats.txt" , emit: sorted_samtools_stats
   
   script:    
   """
   samtools stats ${myBamSample} | grep ^SN | cut -f 2- > ${myBamSample.baseName}_samtools_stats.txt   
   samtools stats ${sorted_mark} | grep ^SN | cut -f 2- > ${sorted_mark.baseName}_samtools_stats.txt
    """
}