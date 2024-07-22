process DNAm_Full_Matrix {

    executor 'local'
    input:
    path files_ch
   
    output:
    val ("DNAm_Full_Matrix.csv"), emit: DNAm_Full_Matrix

    script:
    """
     mkdir -p ${baseDir}/${params.outdir}/DNAm_Full_Matrix
     Rscript ${baseDir}/bin/Full_DNAm_Matrix.R ${params.pipeline} ${baseDir}/${params.outdir}/Methylation ${baseDir}/${params.outdir}/DNAm_Full_Matrix
    """
}

