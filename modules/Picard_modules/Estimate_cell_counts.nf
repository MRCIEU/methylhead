process Estimate_cell_counts {
    
    input:
    
    path full_matrix
    
    publishDir "${params.outdir}/Estimate_cell_count", mode: 'copy' , pattern: '*.csv'
   
    output:
   
    path ("estimate_cell_counts.csv"), optional: true
      
    shell:
    """  
    mkdir -p ${params.outdir}/Estimate_cell_count  
    Rscript ${baseDir}/bin/Estimate_cell_counts.R \
    	    ${params.pipeline} \
	    ${baseDir}/data/blood_cell_types_extended.zip \
	    ${full_matrix} \
	    ${params.outdir}/Estimate_cell_count/estimate_cell_counts.csv
    """
}
