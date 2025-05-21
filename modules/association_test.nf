process association_test {
  
  publishDir "${params.outdir}/association-test", mode: 'copy'
  
  input:
    path qc_files_ch
    path params.phenotype
    path params.models

  output:
    path "association-test-results", emit: association_files

  script:
  """
  Rscript --vanilla ${projectDir}/scripts/test-dnam-assocs.r ${qc_files_ch} ${params.phenotype} ${params.models} association-test-results
  """
}
