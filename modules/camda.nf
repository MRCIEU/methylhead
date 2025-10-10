process camda {
 
  input:
    tuple val(sample_id), path(bam_camda)
    path genome_fasta
    val(samtools_path)
     
  output:   
    tuple val(sample_id), path("${sample_id}_CpG_CAMDA.tsv.gz")     , emit: camda
    tuple val(sample_id), path("${sample_id}_CpG_MethRatio.tsv.gz") , emit: meth
     
   shell:
   """
   SAMTOOLS_PATH=\$(realpath ${samtools_path})
   FASTA=\$(realpath ${genome_fasta})
   python ${projectDir}/scripts/camda.py CAMDA ${bam_camda} \${FASTA} -o ${sample_id} -w ${sample_id} -s \${SAMTOOLS_PATH} -X CG
   gzip ${sample_id}_CpG_CAMDA.tsv
   gzip ${sample_id}_CpG_MethRatio.tsv
   """
}
