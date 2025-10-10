process test_associations {
  
  input:
    path(input_files)
    path(phenotype_file)
    path(model_file)

  output:
    path("associations"), emit: file_path

  script:
  """
  Rscript ${projectDir}/scripts/test-associations.r ${input_files} ${phenotype_file} ${model_file} associations
  """
}
