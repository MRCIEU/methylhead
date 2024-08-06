process Interval_file {

    input:
    path params.panel
    path params.genome_folder

    output:
    path "interval_file", emit: panel

    publishDir "${params.outdir}/Interval_file/", mode: 'copy'

    script:
    """
    picard BedToIntervalList \
    I=${params.panel} \
    O=interval_file \
    SD=${params.genome_folder}.dict
    """
}
