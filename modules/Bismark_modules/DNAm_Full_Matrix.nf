process DNAm_Full_Matrix {

    input:
    path file_ch
    path script

    publishDir "${params.outdir}/DNAm_Full_Matrix/" , mode: 'copy'
    
    output:
    val ("DNAm_Full_Matrix.csv") , emit: DNAm_Full_Matrix
   

    shell:
    """
    Rscript ${script}/scripts/Full_DNAm_Matrix.R ${params.pipeline} ${baseDir}/${params.outdir}/Methylation ${baseDir}/${params.outdir}/DNAm_Full_Matrix
    """ 
}
