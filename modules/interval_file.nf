process interval_file {
  
  publishDir "${params.outdir}/interval-file/", mode: 'copy'

  input:
    path params.target_regions
    path params.genome_folder

  output:
    path "interval_file", emit: target_regions

  script:
    """
    picard BedToIntervalList \
    I=${params.target_regions} \
    O=interval_file \
    SD=${params.genome_folder}
    """
}
