process Methylation_Matrix {

    
    input:
    path files_ch

     publishDir "${params.outdir}/Methylation_Matrix/", mode: 'copy', pattern: '*.pdf, *.csv'

    output:
    val 'Bismark_methylation_matrix.csv' , emit : meth_matrix
    val 'correlation_plot.pdf'    , emit : corr_plot 
    val 'culster_plot.pdf'        , emit : cluster_plot
   
   shell:
    """
    mkdir -p ${baseDir}/${params.outdir}/Methylation_Matrix
    cd ${baseDir}/${params.outdir}/Methylation
    Rscript ${workflow.projectDir}/scripts/Bismark_Methylation_Matrix.R ${baseDir}/${params.outdir}/Methylation_Matrix/
    """
}
