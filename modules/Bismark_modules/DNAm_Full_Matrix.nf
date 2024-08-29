process DNAm_Full_Matrix {

    input:

    path R_files, from: 'R_files/*'
   
    output:
    val "DNAm_Full_Matrix.csv", emit: DNAm_Full_Matrix
      
    publishDir "${params.outdir}/DNAm_Full_Matrix", mode: 'copy'

    shell:
    """
    mkdir -p ${params.outdir}/DNAm_Full_Matrix
    Rscript ${baseDir}/bin/Full_DNAm_Matrix.R ${params.pipeline} ${params.outdir}/Methylation ${params.outdir}/DNAm_Full_Matrix
    """
}



