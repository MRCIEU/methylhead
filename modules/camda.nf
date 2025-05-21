process camda {
 
  publishDir "${params.outdir}/camda-files" , mode: 'copy'
  
  input:
    tuple val(sample_id), path(bam_camda)
     
  output:   
    tuple val(sample_id), path("${sample_id}_CpG_CAMDA.tsv")     , emit: camda
    tuple val(sample_id), path("${sample_id}_CpG_MethRatio.tsv") , emit: meth_ratio
     
   shell:
   """ 
   python ${baseDir}/scripts/camda.py CAMDA ${bam_camda} ${params.genome_folder} -o ${sample_id} -w ${sample_id} -s ${params.samtools_path} -X CG
   """
}
