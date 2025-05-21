process camda_matrix {  
 
  publishDir "${params.outdir}/camda-matrix", mode: 'copy'

  input:
    path(files)
   
  output:
    path ("camda-matrix.csv") , emit: camda_matrix

  shell:
   """
    Rscript --vanilla ${projectDir}/scripts/combine-camda.r ${files} camda-matrix.csv
   """
}
