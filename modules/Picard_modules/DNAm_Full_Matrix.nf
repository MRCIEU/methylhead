process DNAm_Full_Matrix {
    
   executor 'local'
   cpus 16
    
   input:
   path files_ch
   
   output:
   val "DNAm_Full_Matrix.csv", emit: DNAm_Full_Matrix
      
   publishDir "${params.outdir}/DNAm_Full_Matrix", mode: 'copy'

    shell:
    """
    mkdir -p ${params.outdir}/DNAm_Full_Matrix
    Rscript ${baseDir}/bin/Full_DNAm_Matrix.R ${params.pipeline} ${params.outdir}/Methylation ${params.outdir}/DNAm_Full_Matrix
    """
}



