
process Bwa_aligner{

  input:
  tuple val(sample_id), path(fq)
  
  publishDir "${params.outdir}/Bam_files/" , mode: 'copy'
  
  output:
  path("${sample_id}.bam"), emit: bam
      
  script:
  """
  bwameth.py --reference ${params.genome_folder}  ${fq[0]} ${fq[1]} -t 12 | samtools view -b - > ${sample_id}.bam
  """
}
