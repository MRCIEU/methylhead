process Methylation_Matrix {
    executor 'local'   
    
    input:
    path R_files , from: 'R_files/*'         

    publishDir "${params.outdir}/Methylation_Matrix/", mode: 'copy', pattern: '*.pdf, *.csv'
   
    output:
    val  'Methylation_matrix.csv' , emit : meth_matrix
    
    shell:
    """
    mkdir -p ${params.outdir}/Methylation_Matrix
    Rscript ${baseDir}/bin/Methylation_matrix.R ${params.pipeline} ${params.outdir}/Methylation ${params.outdir}/Methylation_Matrix 
    
    """
}
