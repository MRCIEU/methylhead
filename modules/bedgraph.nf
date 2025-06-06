process bedgraph {
 
  publishDir "${params.outdir}/bedgraph-files/", mode: 'copy'
 
  input:    
   tuple val(sample_id), path (sorted_mark)
    
  output:
   tuple val(sample_id), path ("${sample_id}.markdup_CpG.bedGraph") , emit: bedGrap
    
  script:      
     """
     MethylDackel extract --minDepth 10 --maxVariantFrac 0.25 --nOT 0,0,0,98 --nOB 0,0,3,0 --mergeContext ${params.genome_folder} ${sorted_mark} --keepDupes 
     """
}