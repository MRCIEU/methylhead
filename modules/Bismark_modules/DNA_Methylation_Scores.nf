process DNA_Methylation_Scores {

    input:
   
    path full_matrix
     
    publishDir "${params.outdir}/DNA_Methylation_Scores/", mode: 'copy' , pattern: '*.csv'
      
    output:
    path ("DNA_Methylation_Scores.csv"), emit : bismark_scores 
    path ("DNA_Methylation_Sites.csv"), emit : bismark_scores_sites
 
    shell:
    """
    Rscript ${baseDir}/bin/DNA_Methylation_Scores.R ${full_matrix}  DNA_Methylation_Scores.csv DNA_Methylation_Sites.csv
    """
}
