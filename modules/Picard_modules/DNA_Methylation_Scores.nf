
process DNA_Methylation_Scores {

    input:
   
    path Meth_Matrix
     
   output:
    val "Bismark_scores.csv", emit: Bismark_scores
    val "Bismark_smoking_scores.csv", emit: Bismark_smoking_scores

    publishDir "${params.outdir}/DNA_Methylation_Scores/", mode: 'copy'
 
    shell:
    """
    mkdir -p ${params.outdir}/DNA_Methylation_Scores/
    Rscript ${baseDir}/bin/DNA_Methylation_Scores.R ${params.pipeline} ${params.outdir}/Methylation_Matrix ${params.outdir}/DNA_Methylation_Scores
    """
}
