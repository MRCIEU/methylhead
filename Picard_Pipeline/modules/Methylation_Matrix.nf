process Methylation_Matrix {

    input:
    path methylKit

    output:
    val "picard_methylation.csv", emit: meth_matrix


    script:
    """
    mkdir -p ${baseDir}/${params.outdir}/Methylation_Matrix
    cd ${baseDir}/${params.outdir}/methylKit
    Rscript ${workflow.projectDir}/scripts/Picard_Methylation_Matrix.R ${baseDir}/${params.outdir}/Methylation_Matrix/
    """
}