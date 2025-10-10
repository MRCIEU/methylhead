process extract_cpg_methylation {
 
  input:  
    tuple val(sample_id), path (bam_file)
    path  genome_fasta
  
  output:
    tuple val(sample_id), path ("${sample_id}.markdup_CpG.methylKit.gz") , emit: cpg
    
  script:      
  """
  MethylDackel extract --methylKit ${genome_fasta} ${bam_file}
  gzip ${sample_id}.markdup_CpG.methylKit
  """
}
