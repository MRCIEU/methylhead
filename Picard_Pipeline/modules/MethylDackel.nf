
process MethylDackel {
    input:
    path sorted_mark
    
    publishDir "${params.outdir}/MethylDackel/", mode: 'copy'
    
    output:
    path "${sorted_mark.baseName}.markdup_OT.svg" , emit: markdup_OT
    path "${sorted_mark.baseName}.markdup_OB.svg" , emit: markdup_OB
   

   script:    
   """
   MethylDackel mbias ${params.genome_folder} ${sorted_mark} ${sorted_mark.baseName}.markdup --nOT 0,0,0,98 --nOB 0,0,3,0
    """
}