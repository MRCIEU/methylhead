
process Estimate_cell_counts {

    input:
    path full_matrix2

    publishDir "${params.outdir}/Estimate_cell_count/" , mode: 'copy'
    
    output:
    val "Estimate_cell_counts.csv", emit: Estimate_cell_count
   

    shell:
    """
    mkdir -p ${baseDir}/${params.outdir}/Estimate_cell_counts
    Rscript ${workflow.projectDir}/scripts/Estimate_cell_counts.R ${params.pipeline} ${workflow.projectDir}/data ${baseDir}/${params.outdir}/DNAm_Full_Matrix  ${baseDir}/${params.outdir}/Estimate_cell_counts
  
    """
}
