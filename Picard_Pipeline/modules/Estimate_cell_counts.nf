
process Estimate_cell_counts {

    input:
    path data_ch

    publishDir "${params.outdir}/Estimate_cell_count/" , mode: 'copy'
    
    output:
    val "Estimate_cell_counts.csv", emit: Estimate_cell_count
   

    script:
    """
    mkdir -p ${baseDir}/${params.outdir}/Estimate_cell_counts
    cd ${baseDir}/${params.outdir}/methylKit
    Rscript ${workflow.projectDir}/scripts/Estimate_cell_counts.R ${baseDir}/${params.outdir}/Estimate_cell_counts/ ${workflow.projectDir}/data
    """
}