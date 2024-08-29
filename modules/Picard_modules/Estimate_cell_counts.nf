process Estimate_cell_counts {
    
    input:
    
    path R_files, from: 'R_files/*'
    
    publishDir "${params.outdir}/Estimate_cell_count", mode: 'copy' , pattern: '*.csv'
   
    output:
   
    path ("Estimate_cell_counts.csv"), optional: true
      
    shell:
    """  
    mkdir -p ${params.outdir}/Estimate_cell_count
    Rscript ${baseDir}/bin/Estimate_cell_counts.R ${params.pipeline} ${baseDir}/data  ${params.outdir}/Methylation ${params.outdir}/Estimate_cell_count
    """
}
