process dna_methylation_scores {
  
  input:   
    path methylation_matrix
    val assembly
  
  output:
    path ("dna-methylation-scores.csv") , emit : scores
    path ("dna-methylation-sites.csv")  , emit : counts
  
  shell:
    """
    Rscript --vanilla ${projectDir}/scripts/dna-methylation-scores.r ${methylation_matrix}  dna-methylation-scores.csv dna-methylation-sites.csv
    """
}
