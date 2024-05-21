process Samtools_stats {
    input:
    tuple val(sample_id), path(myBamSample)
    tuple val(sample_id), path(sorted_mark)
    
    publishDir "${params.outdir}/Samtools_stats/", mode: 'copy'
    
    output:
    tuple val(sample_id) , path ("${sample_id}_samtools_stats.txt") , emit: mybam_samtools_stats
    tuple val(sample_id) , path ("${sample_id}.markdup_samtools_stats.txt") , emit: sorted_samtools_stats
   
   script:    
   """
   samtools stats ${myBamSample} | grep ^SN | cut -f 2- > ${sample_id}_samtools_stats.txt   
   samtools stats ${sorted_mark} | grep ^SN | cut -f 2- > ${sample_id}.markdup_samtools_stats.txt
    """
}