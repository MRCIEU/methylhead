process DNA_Methylation_Scores {


    input:
   
    path dnascore
     
    publishDir "${params.outdir}/DNA_Methylation_Scores/", mode: 'copy' , pattern: '*.csv'
      
    output:
    path ("Bismark_scores.csv"), emit : Bismark_scores ,optional: true
    path ("Bismark_smoking_scores.csv"), emit : Bismark_smoking_scores, optional: true
 
    shell:
    """
    mkdir -p ${params.outdir}/DNA_Methylation_Scores/
    Rscript ${baseDir}/bin/DNA_Methylation_Scores.R ${params.pipeline} ${params.outdir}/Methylation_Matrix ${params.outdir}/DNA_Methylation_Scores
    """
}
