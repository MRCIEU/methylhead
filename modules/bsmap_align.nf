process bsmap_align {
  
  input:
    tuple val(sample_id), path(fastq)
    path  genome_fasta
    
  output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
      
  script:
  """
  FASTA=\$(realpath ${genome_fasta})
  bsmap -d \${FASTA} -a ${fastq[0]} -b ${fastq[1]} -R -p 16 -o ${sample_id}.bam
  """
}
