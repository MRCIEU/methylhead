
process Methylation_Matrix {

    
     input:
    
     path files_ch
     path script
     
    

     publishDir "${params.outdir}/Methylation_Matrix/", mode: 'copy', pattern: '*.pdf, *.csv'
   
    output:
    val  'Methylation_matrix.csv' , emit : meth_matrix
    val 'correlation_plot.pdf'    , emit : corr_plot 
    val 'culster_plot.pdf'       , emit : cluster_plot
    
    shell:
    """
    Rscript ${script}/scripts/Methylation_matrix.R ${params.pipeline} ${baseDir}/${params.outdir}/Methylation ${baseDir}/${params.outdir}/Methylation_Matrix
    
    """
}
