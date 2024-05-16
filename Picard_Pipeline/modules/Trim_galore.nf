
process Trim_galore {

    input:
    tuple val(sample_id), path(reads)
  
    publishDir "${params.outdir}/Trimmed/" , mode: 'copy'
    
    output:                               
    tuple val(sample_id), path("*_{1,2}.fq.gz") , emit: fq                                  
    tuple val(sample_id), path("*report.txt")                        , emit: log     , optional: true 
 
    script:
    def trim_galore_cmd = "trim_galore --paired ${reads[0]} ${reads[1]} --gzip"
    "${trim_galore_cmd}"

}   