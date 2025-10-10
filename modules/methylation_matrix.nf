process methylation_matrix {
  
  input:
    path(files)
    val(assembly)
  
  output:
    path ("methylation-matrix.csv.gz") , emit: meth
    path ("coverage-matrix.csv.gz")    , emit: coverage
   
  shell:
    """
    Rscript ${projectDir}/scripts/methylation-matrix.r ${files} methylation-matrix.csv coverage-matrix.csv ${assembly}
    gzip methylation-matrix.csv
    gzip coverage-matrix.csv
    """
}


