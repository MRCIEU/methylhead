
process DNA_Methylation_Scores {

    input:
   
    path Meth_Matrix

    publishDir "${params.outdir}/Picard_scores/" , mode: 'copy'
    
    output:
    val "Picard_scores.csv", emit: Picard_scores
    val "Picard_smoking_scores.csv" , emit : Picard_smoking_scores
   

    shell:
    """
    mkdir -p ${baseDir}/${params.outdir}/Picard_scores 
    cd ${baseDir}/${params.outdir}/Methylation_Matrix/
    Rscript ${workflow.projectDir}/scripts/Picard_scores.R ${baseDir}/${params.outdir}/Picard_scores/
    """
}