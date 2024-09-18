process Illumina_Matrix {
      
    
    input:
    
    path full_matrix         

    publishDir "${params.outdir}/Methylation_Matrix/" , mode: 'copy' , pattern: '*.csv'
   

    output:
    path 'Methylation_matrix.csv', emit: meth_matrix    
    
    shell:   
    """
    Rscript ${baseDir}/bin/Methylation_matrix.R ${full_matrix} Illumina_matrix.csv
    """
}



