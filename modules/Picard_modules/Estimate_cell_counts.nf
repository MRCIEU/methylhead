process Estimate_cell_counts {

         
   executor 'local'
   cpus 16


    input:
    path full_matrix2
   
    output:
   
    val "Estimate_cell_counts.csv", emit: Estimate_cell_count
   
    publishDir "${params.outdir}/Estimate_cell_count", mode: 'copy'
    
    shell:
    """  
    mkdir -p ${params.outdir}/Estimate_cell_count
    Rscript ${baseDir}/bin/Estimate_cell_counts.R ${params.pipeline} ${baseDir}/data ${params.outdir}/DNAm_Full_Matrix ${params.outdir}/Estimate_cell_count
     """
}
