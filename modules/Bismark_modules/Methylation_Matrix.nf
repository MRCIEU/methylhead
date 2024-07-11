process Methylation_Matrix {

     
     
     input:
   
     path R_files , from: 'R_files/*'         

     publishDir "${params.outdir}/Methylation_Matrix/", mode: 'copy', pattern: '*.pdf, *.csv'
   
    output:
    val  'Methylation_matrix.csv' , emit : meth_matrix
    val 'correlation_plot.pdf'    , emit : corr_plot 
    val 'culster_plot.pdf'       , emit : cluster_plot
    
    shell:
    """
    mkdir -p ${baseDir}/${params.outdir}/Methylation_Matrix
    Rscript ${baseDir}/bin/Methylation_matrix.R ${params.pipeline} ${baseDir}/${params.outdir}/Methylation ${baseDir}/${params.outdir}/Methylation_Matrix 
    
    """
}
