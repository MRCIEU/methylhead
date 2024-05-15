
process Estimate_cell_counts {

    input:
    path full_matrix2

    publishDir "${params.outdir}/Estimate_cell_count/" , mode: 'copy'
    
    output:
    val "Estimate_cell_counts.csv", emit: Estimate_cell_count
   

    shell:
    """
    mkdir -p ${baseDir}/${params.outdir}/Estimate_cell_counts
    cd ${baseDir}/${params.outdir}/DNAm_Full_Matrix
    Rscript ${workflow.projectDir}/scripts/Estimate_cell_counts.R ${baseDir}/${params.outdir}/Estimate_cell_counts/ ${workflow.projectDir}/data
    """
}
