process DNAm_Full_Matrix {

    input:
    path file_ch

    publishDir "${params.outdir}/DNAm_Full_Matrix/" , mode: 'copy'
    
    output:
    val ("DNAm_Full_Matrix.csv") , emit: DNAm_Full_Matrix
   

   shell:
    """
    mkdir -p ${baseDir}/${params.outdir}/DNAm_Full_Matrix
    cd ${baseDir}/${params.outdir}/Methylation
    Rscript ${workflow.projectDir}/scripts/Bismark_full_DNAm_matrix.R ${baseDir}/${params.outdir}/DNAm_Full_Matrix/ 
    """

}