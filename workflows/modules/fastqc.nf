
process fastqc {
  
  publishDir "${params.outdir}/fastqc-files/" , mode: 'copy'
   
  input:
    tuple val(sample_id), path(reads)
 
  output:
    tuple val(sample_id), path("*.html") , emit: html
    tuple val(sample_id), path("*.zip")  , emit: zip
    
  script:
    def fastqc_cmd = "fastqc ${reads[0]} ${reads[1]}"  
    "${fastqc_cmd}"
}
