process DNAm_Matrix {
   
    input:
    path(files)

    publishDir "${params.outdir}/DNAm_Matrix", mode: 'copy'

    output:

    path ("Methylation_matrix.csv") , emit: meth_matrix
    path ("coverage_matrix.csv")    , emit: coverage_matrix
   
    shell:
    """
    Rscript ${baseDir}/bin/Full_DNAm_Matrix.R ${params.pipeline} ${files} Methylation_matrix.csv coverage_matrix.csv
    """
}
