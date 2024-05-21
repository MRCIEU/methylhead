process MethylKit{

    input:  
    tuple val(sample_id), path (sorted_mark)
    
    publishDir "${params.outdir}/Methylation/", mode: 'copy'
    
    output:
     tuple val(sample_id), path ("${sample_id}.markdup_CpG.methylKit") , emit: methylKit
    script:      
     """
    MethylDackel extract --methylKit ${params.genome_folder} ${sorted_mark} 
     """
}
