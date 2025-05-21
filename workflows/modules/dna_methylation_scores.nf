process dna_methylation_scores {
  
  publishDir "${params.outdir}/dna-methylation-scores/", mode: 'copy' , pattern: '*.csv'
    
  input:   
    path methylation_matrix
    
  output:
    path ("dna-methylation-scores.csv") , emit : dna_methylation_scores
    path ("dna-methylation-sites.csv")  , emit : sites_scores_sites
  
  shell:
    """
    Rscript --vanilla ${projectDir}/scripts/dna-methylation-scores.r ${methylation_matrix}  dna-methylation-scores.csv dna-methylation-sites.csv
    """
}
