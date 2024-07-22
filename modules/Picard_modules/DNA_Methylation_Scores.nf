
process DNA_Methylation_Scores {
    executor 'local'
    input:   
    path Meth_Matrix

    publishDir "${params.outdir}/DNA_Methylation_Scores/" , mode: 'copy'
    
    output:
    val ("Bismark_scores.csv"), emit: Bismark_scores
    val ("Bismark_smoking_scores.csv") , emit : Bismark_smoking_scores
   

    shell:
    """
    mkdir -p ${baseDir}/${params.outdir}/DNA_Methylation_Scores/
    Rscript ${baseDir}/bin/DNA_Methylation_Scores.R ${params.pipeline} ${baseDir}/${params.outdir}/Methylation_Matrix ${baseDir}/${params.outdir}/DNA_Methylation_Scores
    """
}
