process Illumina_Matrix {
          
    input:

    path full_matrix         
    publishDir "${params.outdir}/Illumina_Matrix/" , mode: 'copy' , pattern: '*.csv'
   

    output:
    path 'Illumina_matrix.csv', emit: Illumina_matrix    
    
    shell:   
    """
    Rscript ${baseDir}/bin/Methylation_matrix.R ${full_matrix} Illumina_matrix.csv
    """
}



