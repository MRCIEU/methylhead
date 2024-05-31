process BSmap_Aligment {

      input:
      tuple val(sample_id), path(fq)
  
      publishDir "${params.outdir}/BSmap_Aligment/" , mode: 'copy'
  
      output:
      tuple val(sample_id), path("${sample_id}.bam"), emit: bam
      
      script:
      """
      bsmap -d ${params.genome_folder} -a ${fq[0]} -b ${fq[1]} -R -p 16 -o ${sample_id}.bam
      """

}
