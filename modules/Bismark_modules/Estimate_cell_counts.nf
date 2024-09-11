process Estimate_cell_counts {
    
    input:
    
    path full_matrix
    
    publishDir "${params.outdir}/Estimate_cell_count", mode: 'copy' , pattern: '*.csv'
   
    output:
   
    path ("estimate_cell_counts.csv") , emit: estimate_cell_counts
      
    shell:
    """   
    Rscript ${baseDir}/bin/Estimate_cell_counts.R ${baseDir}/data/blood_cell_types_extended.zip ${full_matrix} estimate_cell_counts.csv
    """
}
