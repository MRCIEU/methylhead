process DNAm_Full_Matrix {
   
    input:
    path(files)

    publishDir "${params.outdir}/DNAm_Full_Matrix", mode: 'copy'

    output:

    path ("Methylation_matrix.csv") , emit: meth_matrix
    path ("coverage_matrix.csv"), emit: covarage_matrix
   
    shell:
    """
    Rscript ${baseDir}/bin/Full_DNAm_Matrix.R ${params.pipeline} ${files} Methylation_matrix.csv coverage_matrix.csv
    """
}
