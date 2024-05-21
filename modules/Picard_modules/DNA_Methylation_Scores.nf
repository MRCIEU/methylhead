
process DNA_Methylation_Scores {

    input:
   
    path Meth_Matrix

    publishDir "${params.outdir}/DNA_Methylation_Scores/" , mode: 'copy'
    
    output:
    val "Picard_scores.csv", emit: Picard_scores
    val "Picard_smoking_scores.csv" , emit : Picard_smoking_scores
   

   
    shell:
    """
   mkdir -p ${baseDir}/${params.outdir}/DNA_Methylation_Scores/
    Rscript ${workflow.projectDir}/scripts/DNA_Methylation_Scores.R ${params.pipeline} ${baseDir}/${params.outdir}/Methylation_Matrix  ${baseDir}/${params.outdir}/DNA_Methylation_Scores
    """
}