process interval_file {
  
  publishDir "${params.outdir}/interval-file/", mode: 'copy'

  input:
    path params.panel
    path params.genome_folder

  output:
    path "interval_file", emit: panel

  script:
    """
    picard BedToIntervalList \
    I=${params.panel} \
    O=interval_file \
    SD=${params.genome_folder}.dict
    """
}
