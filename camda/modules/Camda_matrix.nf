process Camda_matrix {
   
    input:
    path(files)
   

    publishDir "${params.outdir}/Camda_matrix", mode: 'copy'

    output:
   
    path ("Camda_matrix.csv") , emit: camda_matrix
   
    shell:
    """
    Rscript ${workflow.projectDir}/scripts/combinecamda.R ${files} Camda_matrix.csv
    """
}
