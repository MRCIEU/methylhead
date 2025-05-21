process methylation_matrix_process {
  
  publishDir "${params.outdir}/methylation-matrix", mode: 'copy'  
  
  input:
    path(files)
  
  output:
    path ("methylation-matrix.csv") , emit: meth_matrix
    path ("coverage-matrix.csv")    , emit: coverage_matrix
   
  shell:
    """
    Rscript ${baseDir}/scripts/methylation-matrix.r ${files} methylation-matrix.csv coverage-matrix.csv
    """
}
