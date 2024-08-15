process Methylation_Matrix {
      
    
    input:
    
    path R_files , from: 'R_files/*'         

    publishDir "${params.outdir}/Methylation_Matrix/", mode: 'move' , pattern: '*.csv' , saveAs: {  Methylation_matrix.csv } 
   

    output:
    path 'Methylation_matrix.csv', emit: meth_matrix  , optional: true   
    
    shell:   
    """
    mkdir -p ${params.outdir}/Methylation_Matrix
    Rscript ${baseDir}/bin/Methylation_matrix.R ${params.pipeline} ${params.outdir}/Methylation  ${params.outdir}/Methylation_Matrix
    """
}


