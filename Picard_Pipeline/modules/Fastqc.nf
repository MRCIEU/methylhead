
process Fastqc {
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip") , emit: zip
    
    publishDir "${params.outdir}/fastqc/" , mode: 'copy'
    
    script:
    def fastqc_cmd = "fastqc ${reads[0]} ${reads[1]}"  
    "${fastqc_cmd}"

}