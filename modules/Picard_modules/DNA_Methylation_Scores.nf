process DNA_Methylation_Scores {

    input:
   
    path R_files , from: 'R_files/*'
     
    publishDir "${params.outdir}/DNA_Methylation_Scores/", mode: 'copy' , pattern: '*.csv'
      
    output:
    path ("DNA_Methylation_Scores.csv"), emit : Picard_scores ,optional: true
    path ("Picard_smoking_scores.csv"), emit : Picard_smoking_scores, optional: true
 
    shell:
    """
    mkdir -p ${params.outdir}/DNA_Methylation_Scores/
    Rscript ${baseDir}/bin/DNA_Methylation_Scores.R ${params.pipeline} ${params.outdir}/Methylation ${params.outdir}/DNA_Methylation_Scores
    """
}
