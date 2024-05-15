
process DNA_Methylation_Scores {

    input:
   
    path Meth_Matrix

    publishDir "${params.outdir}/DNA_Methylation_Scores/" , mode: 'copy'
    
    output:
    val ("Bismark_scores.csv"), emit: Bismark_scores
    val ("Bismark_smoking_scores.csv") , emit : Bismark_smoking_scores
   

    shell:
    """
    mkdir -p ${baseDir}/${params.outdir}/DNA_Methylation_Scores/
    cd ${baseDir}/${params.outdir}/Methylation_Matrix/
    Rscript ${workflow.projectDir}/scripts/Bismark_scores.R ${baseDir}/${params.outdir}/DNA_Methylation_Scores/
    """
}