process EWAS_Analysis {  
    
    container "/home/ag24712/pipeline_update/meth_analysis.sif"

    input:
    path(files)
    path(phenotype)
    path(models)
   
    publishDir "${params.outdir}/EWAS_Analysis", mode: 'copy'

    output:
    path "output/**", emit: ewas_results 

    shell:
    """
    Rscript ${baseDir}/bin/test-dnam-assocs.r ${params.outdir} ${params.phenotype} ${params.models} ${params.removed_samples} output
    """
}



