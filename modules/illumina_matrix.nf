process illumina_matrix {
  
  input:
   path methylation_matrix 
   val assembly        

  output:
   path 'illumina-matrix.csv.gz', emit: meth
    
  shell:   
   """
   Rscript --vanilla ${projectDir}/scripts/illumina-matrix.r ${methylation_matrix} illumina-matrix.csv $assembly
   gzip illumina-matrix.csv
   """
}



