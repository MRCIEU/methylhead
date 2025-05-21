process illumina_matrix_450k {
  
  publishDir "${params.outdir}/illumina-matrix-450k/" , mode: 'copy' , pattern: '*.csv'

  input:
   path methylation_matrix         
    
  output:
   path 'illumina-matrix.csv', emit: illumina_matrix    
    
  shell:   
   """
   Rscript --vanilla ${projectDir}/scripts/illumina-matrix.r ${methylation_matrix} illumina-matrix.csv
   """
}



