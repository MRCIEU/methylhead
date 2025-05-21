process estimate_cell_counts {
  
  publishDir "${params.outdir}/estimate-cell-counts", mode: 'copy' , pattern: '*.csv'   
  
  input: 
    path methylation_matrix
    
  output:
    path ("estimate-cell-counts.csv") , emit: estimate_cell_counts
      
  shell:
    """   
    Rscript --vanilla ${projectDir}/scripts/estimate-cell-counts.r ${baseDir}/data/blood_cell_types_extended.zip ${methylation_matrix} estimate-cell-counts.csv
    """
}
