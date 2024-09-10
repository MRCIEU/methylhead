process DNAm_Full_Matrix {

    input:
    
    path sample_meth_files, from: 'sample_meth_files/*'

    publishDir "${params.outdir}/DNAm_Full_Matrix", mode: 'copy', pattern: '*.csv' 

    output: 
    
    path ("DNAm_Full_Matrix.csv") , emit: DNAm_Full_Matrix
    
    shell:
    """
    Rscript ${baseDir}/bin/Full_DNAm_Matrix.R ${params.pipeline} ${sample_meth_files} ./
    """
}

