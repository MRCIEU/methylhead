
process Estimate_cell_counts {

    input:
    path R_files, from: 'R_files/*'
   
    output:
   
    val "Estimate_cell_counts.csv", emit: Estimate_cell_count
   
    publishDir "${params.outdir}/Estimate_cell_count", mode: 'copy'
    
    shell:
    """  
    mkdir -p ${params.outdir}/Estimate_cell_count
    Rscript ${baseDir}/bin/Estimate_cell_counts.R ${params.pipeline} ${baseDir}/data ${params.outdir}/Methylation_CpG_only ${params.outdir}/Estimate_cell_count
     """
}
