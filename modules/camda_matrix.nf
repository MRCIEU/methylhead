process camda_matrix {  
 
  input:
    path(files)
    val assembly
   
  output:
    path ("camda-matrix.csv.gz") , emit: scores

  shell:
   """
    Rscript --vanilla ${projectDir}/scripts/camda-matrix.r ${files} camda-matrix.csv $assembly
    gzip camda-matrix.csv
   """
}
