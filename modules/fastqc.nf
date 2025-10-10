process fastqc {
  
  input:
    tuple val(sample_id), path(reads)
 
  output:
    tuple val(sample_id), path("*.html") , emit: html
    tuple val(sample_id), path("*.zip")  , emit: zip
    
  script:
  """      
  fastqc ${reads[0]} ${reads[1]}
  """
}
