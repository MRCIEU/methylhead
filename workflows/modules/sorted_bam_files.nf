process sorted_bam_files {

  publishDir "${params.outdir}/sorted-bam-files/" , mode: 'copy'

  input:
      tuple val(sample_id),  path (sortedbam)

  output:   
      tuple val(sample_id) , path ("${sample_id}_sorted.bam")     , emit: sorted_bam
      tuple val(sample_id) , path ("${sample_id}_sorted.bam.bai") , emit: bai
     
  script:
    """
     samtools sort -t 16 -m 16GiB -l 0 ${sortedbam} -o ${sample_id}_sorted.bam     
     samtools index -@ 16 ${sample_id}_sorted.bam  > ${sample_id}_sorted.bam.bai
    """
}
