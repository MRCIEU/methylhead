process Association_test {
  
  publishDir "${params.outdir}/Association_test", mode: 'copy'
  
  input:
    path qc_files_ch
    path params.phenotype
    path params.models

  output:
    path "Association_test_results", emit: association_files

  script:
  """
  Rscript ${baseDir}/bin/test-dnam-assocs.r ${qc_files_ch} ${params.phenotype} ${params.models} Association_test_results
  """
}
