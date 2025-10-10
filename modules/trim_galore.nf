process trim_galore {

  input:
      tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("*_{1,2}.fq.gz") , emit: fq
    tuple val(sample_id), path("*report.txt")   , emit: log   , optional: true
    tuple val(sample_id), path("*.html")        , emit: html  , optional: true  

  script:
  """
  trim_galore --paired ${reads[0]} ${reads[1]} --gzip
  """
}


