process bsmap_aligment {
  publishDir "${params.outdir}/bsmap-aligment-bam-files/" , mode: 'copy'
  
  input:
    tuple val(sample_id), path(fq)
    path  genome_fasta
    
  output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
      
  script:
   """
   bsmap -d ${params.genome_fasta} -a ${fq[0]} -b ${fq[1]} -R -p 16 -o ${sample_id}.bam
   """
}
