process Methylation_Matrix {

    
    input:
    path files_ch

    output:
    val "Bismark_methylation.csv" , emit : meth_matrix
   
    publishDir "${params.outdir}/Methylation_Matrix/" , mode: 'copy'

    script:
    """
    mkdir -p ${baseDir}/${params.outdir}/Methylation_Matrix
    cd ${baseDir}/${params.outdir}/Methylation
    Rscript ${workflow.projectDir}/scripts/Bismark_Methylation_Matrix.R ${baseDir}/${params.outdir}/Methylation_Matrix/
    """
}
