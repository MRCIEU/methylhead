process samtools_stats {
  
  publishDir "${params.outdir}/samtools-stats/"    , mode: 'copy'

  input:
    tuple val(sample_id), path(mybamsample)
    tuple val(sample_id), path(sorted_mark)
      
  output:
    tuple val(sample_id) , path ("${sample_id}_samtools_stats.txt")         , emit: mybam_samtools_stats
    tuple val(sample_id) , path ("${sample_id}.markdup_samtools_stats.txt") , emit: sorted_samtools_stats
   
  script:    
   """
   samtools stats ${mybamsample} | grep ^SN | cut -f 2- > ${sample_id}_samtools_stats.txt   
   samtools stats ${sorted_mark} | grep ^SN | cut -f 2- > ${sample_id}.markdup_samtools_stats.txt
    """
}