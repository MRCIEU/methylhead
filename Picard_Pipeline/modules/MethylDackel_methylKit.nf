process MethylDackel_methylKit{

    input:  
    path sorted_mark
    
    publishDir "${params.outdir}/methylKit/", mode: 'copy'
    
    output:
    path "${sorted_mark.baseName}_CpG.methylKit" , emit: methylKit
    
    script:      
     """
    MethylDackel extract --methylKit ${params.genome_folder} ${sorted_mark} 
     """
}
