process DNAm_Full_Matrix {

    input:
    path R_files, from: 'R_files/*'

    publishDir "${params.outdir}/DNAm_Full_Matrix", mode: 'copy', pattern: '*.csv'

    output: 
    path("DNAm_Full_Matrix.csv"), optional: true , emit: DNAm_Full_Matrix
    
    shell:
    """
    mkdir -p ${params.outdir}/DNAm_Full_Matrix
    Rscript ${baseDir}/bin/Full_DNAm_Matrix.R ${params.pipeline} ${params.outdir}/Methylation ${params.outdir}/DNAm_Full_Matrix
    """
}
