process sort_bam {

  input:
      tuple val(sample_id),  path (bamfile)

  output:   
      tuple val(sample_id) , path ("${sample_id}_sorted.bam")     , emit: bam
      tuple val(sample_id) , path ("${sample_id}_sorted.bam.bai") , emit: bai
     
  script:
    """
     samtools sort -t 16 -m 16GiB -l 0 ${bamfile} -o ${sample_id}_sorted.bam     
     samtools index -@ 16 ${sample_id}_sorted.bam  > ${sample_id}_sorted.bam.bai
    """
}
