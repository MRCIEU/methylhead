process DNAm_Full_Matrix {
    input:
    path R_files , from: 'R_files/*'
   
    output:
    val ("DNAm_Full_Matrix.csv"), emit: DNAm_Full_Matrix

    script:
    """
     mkdir -p ${baseDir}/${params.outdir}/DNAm_Full_Matrix
     Rscript ${baseDir}/bin/Full_DNAm_Matrix.R ${params.pipeline} ${baseDir}/${params.outdir}/Methylation ${baseDir}/${params.outdir}/DNAm_Full_Matrix
    """
}

