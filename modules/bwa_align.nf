process bwa_align {
  
  input:
   tuple val(sample_id), path(fq)
   path  genome_fasta
  
  output:
   tuple val(sample_id), path("${sample_id}.bam"), emit: bam
      
  script:
  """
  FASTA=\$(realpath ${genome_fasta})
  bwameth.py --reference \${FASTA} ${fq[0]} ${fq[1]} -t 12 | samtools view -b - > ${sample_id}.bam
  """
}
