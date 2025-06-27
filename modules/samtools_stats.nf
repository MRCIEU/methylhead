process samtools_stats {

    publishDir "${params.outdir}/samtools-stats/", mode: 'copy'

    input:
        tuple val(sample_id), path(mybamsample)
        tuple val(sample_id), path(sorted_mark)
        tuple val(sample_id), path(sorted_ch)
        tuple val(sample_id), path(sorted_ch_bai)
        path panel

    output:
        tuple val(sample_id), path("${sample_id}_samtools_stats.txt")        , emit: mybam_samtools_stats
        tuple val(sample_id), path("${sample_id}_markdup_samtools_stats.txt"), emit: sorted_samtools_stats
        tuple val(sample_id), path("${sample_id}_panel_depth.tsv")           , emit: panel_depth
        tuple val(sample_id), path("${sample_id}_read_counts.tsv")           , emit: read_counts

    script:
    """
    # ---------- basic QC ----------
    samtools stats "${mybamsample}" | awk '/^SN/ {print \$2, \$3}' > "${sample_id}_samtools_stats.txt"
    samtools stats "${sorted_mark}" | awk '/^SN/ {print \$2, \$3}' > "${sample_id}_markdup_samtools_stats.txt"

    # ---------- panel depth ----------
    total_depth=\$(samtools bedcov "${panel}" "${sorted_ch}" | \
              awk '{sum+=\$4} END {print sum}')

    printf "Sample_id\tTotal_Depth\n%s\t%s\n" \
       "${sample_id}" "\${total_depth}" \
       > "${sample_id}_panel_depth.tsv"
    
    # ---------- read counts ----------
    
    total_reads=\$(samtools view -c -F 4 -q 20 -F 1024  "${sorted_ch}")
    panel_reads=\$(samtools view -c -F 4 -q 20 -F 1024 -L "${panel}" "${sorted_ch}")

    ## define a Bash variable from the channel value
    sample_id_bash="${sample_id}"

    ## strip the literal suffix '_sorted'
    clean_id="\${sample_id_bash%_sorted}"

    printf "SampleID\tTotal_mapped\tPanel_mapped\n%s\t%s\t%s\n" \
       "\${clean_id}" "\${total_reads}" "\${panel_reads}" \
       > "${sample_id}_read_counts.tsv"
    """
}
