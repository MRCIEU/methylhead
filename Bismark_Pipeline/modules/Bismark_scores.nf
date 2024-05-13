
process Bismark_scores {

    input:
   
    path Meth_Matrix

    publishDir "${params.outdir}/Bismark_scores/" , mode: 'copy'
    
    output:
    val "Bismark_scores.csv", emit: Bismark_scores
    val "Bismark_smoking_scores.csv" , emit : Bismark_smoking_scores
   

    script:
    """
    mkdir -p ${baseDir}/${params.outdir}/Bismark_scores/
    cd ${baseDir}/${params.outdir}/Methylation_Matrix/
    Rscript ${workflow.projectDir}/scripts/Bismark_scores.R ${baseDir}/${params.outdir}/Bismark_scores/
    """
}