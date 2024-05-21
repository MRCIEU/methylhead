process Methylation_Matrix {

    input:
    path methylKit

    output:
    val  'Methylation_matrix.csv' , emit : meth_matrix
    val 'correlation_plot.pdf'    , emit : corr_plot 
    val 'culster_plot.pdf'       , emit : cluster_plot
    
    
    publishDir "${params.outdir}/Methylation_Matrix/" , mode: 'copy'

   

    shell:
    """
    mkdir -p ${baseDir}/${params.outdir}/Methylation_Matrix
    Rscript ${workflow.projectDir}/scripts/Methylation_matrix.R ${params.pipeline} ${baseDir}/${params.outdir}/Methylation ${baseDir}/${params.outdir}/Methylation_Matrix
    """
}
