process methylkit{
 
  publishDir "${params.outdir}/methylation-files/", mode: 'copy'

  input:  
    tuple val(sample_id), path (sorted_mark)
    path  genome_folder
  
  output:
    tuple val(sample_id), path ("${sample_id}.markdup_CpG.methylKit") , emit: methylKit_CpG
    tuple val(sample_id), path ("${sample_id}.markdup_CHH.methylKit") , emit: methylKit_CHH
    tuple val(sample_id), path ("${sample_id}.markdup_CHG.methylKit") , emit: methylKit_CHG   
    
  script:      
  """
  MethylDackel extract --CHH --CHG --methylKit ${params.genome_folder} ${sorted_mark} 
  """
}
