process extract_cpg_bedgraph {
 
  input:    
   tuple val(sample_id), path (sorted_mark)
   path  genome_fasta
    
  output:
   tuple val(sample_id), path ("${sample_id}.markdup_CpG.bedGraph.gz") , emit: bedgraph
    
  script:      
     """
     MethylDackel extract --minDepth 10 --maxVariantFrac 0.25 --nOT 0,0,0,98 --nOB 0,0,3,0 --mergeContext ${genome_fasta} ${sorted_mark} --keepDupes
     gzip ${sample_id}.markdup_CpG.bedGraph
     """
}