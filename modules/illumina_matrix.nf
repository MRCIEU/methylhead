process illumina_matrix {
  
  publishDir "${params.outdir}/illumina-matrix/" , mode: 'copy' , pattern: '*.csv'

  input:
   path methylation_matrix 
   val assembly        

  output:
   path 'illumina-matrix.csv', emit: illumina_matrix    
    
  shell:   
   """
   Rscript --vanilla ${projectDir}/scripts/illumina-matrix.r ${methylation_matrix} illumina-matrix.csv $assembly
   """
}



