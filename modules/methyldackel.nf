process methyldackel {
  
  publishDir "${params.outdir}/methyldackel-files/", mode: 'copy'
  
  input:
    tuple val(sample_id), path (sorted_mark)
    path reference 
    
  output:
    tuple val(sample_id), path ("${sample_id}.markdup_OT.svg") , emit: markdup_OT
    tuple val(sample_id), path ("${sample_id}.markdup_OB.svg") , emit: markdup_OB
   
  script:    
  """
  MethylDackel mbias ${reference} ${sorted_mark} ${sample_id}.markdup --nOT 0,0,0,98 --nOB 0,0,3,0
  """
}
