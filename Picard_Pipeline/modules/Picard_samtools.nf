
process Picard_samtools {

    input:
    path sorted_mark

    publishDir "${params.outdir}/picard/", mode: 'copy'

    output:
    path "${sorted_mark}.bai", emit: sorted_bai

    script:
    """
    samtools index -@ 16 ${sorted_mark}
    """
}
