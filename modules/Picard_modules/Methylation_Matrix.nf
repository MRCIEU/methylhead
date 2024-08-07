process Methylation_Matrix {

     
     
     input:
   
     path files_ch         

        
    output:
    
    val  "Methylation_matrix.csv"  , emit : meth_matrix
    val  "correlation_plot.pdf"    , emit : corr_plot 
    val  "culster_plot.pdf"        , emit : cluster_plot
    
    publishDir "${params.outdir}/Methylation_Matrix/", mode: 'copy', pattern: '*.pdf, *.csv'
    
    shell:
    """
    mkdir -p ${params.outdir}/Methylation_Matrix
    Rscript ${baseDir}/bin/Methylation_matrix.R ${params.pipeline} ${params.outdir}/Methylation ${params.outdir}/Methylation_Matrix 
    
    """
}
