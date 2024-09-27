process Camda_matrix {  
    input:
    path(files)
   
    publishDir "${params.outdir}/camda_matrix", mode: 'copy'

    output:
    path ("camda_matrix.csv") , emit: camda_matrix

    shell:
    """
    Rscript ${baseDir}/bin/combinecamda.R ${files} camda_matrix.csv
    """
}
