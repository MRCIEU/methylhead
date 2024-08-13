process Methylation_Matrix {

     cpus 16
     
     input:
   
     path files_ch         

        
    output:
    
    val  "Methylation_matrix.csv"  , emit : meth_matrix

    publishDir "${params.outdir}/Methylation_Matrix/", mode: 'copy', pattern: '*.pdf, *.csv'
    
    shell:
    """
    mkdir -p ${params.outdir}/Methylation_Matrix
    Rscript ${baseDir}/bin/Methylation_matrix.R ${params.pipeline} ${params.outdir}/Methylation ${params.outdir}/Methylation_Matrix 
    
    """
}
