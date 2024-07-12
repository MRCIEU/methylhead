
process Alignment {

 input:
  tuple val(sample_id), path(fq)
  path  genome_folder

  publishDir "${params.outdir}/Bam_files/" , mode: 'copy'
  
  output:
  tuple val(sample_id), path("${sample_id}.bam"), emit: bam
      
 script:
    def bwa_meth_cmd = "bwameth.py --reference ${params.genome_folder} ${fq[0]} ${fq[1]} -t 12 | samtools view -b - > ${sample_id}.bam"
    "${bwa_meth_cmd}"		
}
