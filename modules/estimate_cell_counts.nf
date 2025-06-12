process estimate_cell_counts {
  
  publishDir "${params.outdir}/estimate-cell-counts", mode: 'copy' , pattern: '*.csv'   
  
  input: 
    path methylation_matrix
    path cell_reference
    
  output:
    path ("estimate-cell-counts.csv") , emit: estimate_cell_counts
      
  shell:
    """   
    Rscript --vanilla ${projectDir}/scripts/estimate-cell-counts.r ${cell_reference} ${methylation_matrix} estimate-cell-counts.csv
    """
}
