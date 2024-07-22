
process Fastqc {
    executor 'local'  
    
    input:
    tuple val(sample_id), path(reads)
    val(t_param)
    val(memory_param)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip") , emit: zip
    
    publishDir "${params.outdir}/Fastqc/" , mode: 'copy'
    
    script:
    def fastqc_cmd = "fastqc ${reads[0]} ${reads[1]} --output_dir ${params.outdir}"
    if (t_param) {
        fastqc_cmd += " --threads ${params.t_param}"
    }
    if (memory_param) {
        fastqc_cmd += " --memory ${params.memory_param}"
    }    
    "${fastqc_cmd}"


}
