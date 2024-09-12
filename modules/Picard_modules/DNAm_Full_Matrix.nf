process DNAm_Full_Matrix {
   
    input:
    path(files)
   

    publishDir "${params.outdir}/DNAm_Full_Matrix", mode: 'copy'

    output:
   
    path ("DNAm_Full_Matrix.csv") , emit: DNAm_Full_Matrix
   
    shell:
    """
    Rscript ${baseDir}/bin/Full_DNAm_Matrix.R ${params.pipeline} ${files} DNAm_Full_Matrix.csv
    """
}
