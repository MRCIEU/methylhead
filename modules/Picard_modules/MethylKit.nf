process MethylKit{

    input:  
    tuple val(sample_id), path (sorted_mark)
    
    publishDir "${params.outdir}/Methylation/", mode: 'copy'
    
    output:
     tuple val(sample_id), path ("${sample_id}.markdup_CpG.methylKit") , emit: methylKit
     tuple val(sample_id), path ("${sample_id}.markdup_CHH.methylKit") , emit: methylKit_CHH
     tuple val(sample_id), path ("${sample_id}.markdup_CHG.methylKit") , emit: methylKit_CHG   
    
    script:      
     """
    MethylDackel extract --CHH --CHG --methylKit ${params.genome_folder} ${sorted_mark} 
     """
}
