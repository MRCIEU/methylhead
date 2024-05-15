process Trim_galore {
  input:
    tuple val(sample_id), path(reads)
    val(cores)
    
    output:                               
    tuple val(sample_id), path("*_{1,2}.fq.gz") , emit: fq                                  
    tuple val(sample_id), path("*report.txt")                        , emit: log     , optional: true
    tuple val(sample_id), path("*.html")                             , emit: html    , optional: true   
   
    publishDir "${params.outdir}/Trimmed/" , mode: 'copy'
    
 
    script:
    def trim_galore_cmd = "trim_galore --paired ${reads[0]} ${reads[1]} --gzip"
    if(cores) {
       trim_galore_cmd += " --cores ${cores}"
    }
    "${trim_galore_cmd}"

}  