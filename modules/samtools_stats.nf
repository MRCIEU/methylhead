process samtools_stats {

    publishDir "${params.outdir}/samtools-stats/", mode: 'copy'

    input:
        tuple val(sample_id), path(mybamsample)
        tuple val(sample_id), path(sorted_mark)
        tuple val(sample_id), path(sorted_ch)
        tuple val(sample_id), path(sorted_ch_bai)
        path panel

    output:
        tuple val(sample_id), path("${sample_id}-samtools-stats.txt")        , emit: mybam_samtools_stats
        tuple val(sample_id), path("${sample_id}-markdup-samtools-stats.txt"), emit: sorted_samtools_stats
        tuple val(sample_id), path("${sample_id}-read-counts.tsv")           , emit: read_counts

    script:
    """
    # ---------- basic QC ----------
    samtools stats "${mybamsample}" | awk '/^SN/ {print \$2, \$3}' > "${sample_id}-samtools-stats.txt"
    samtools stats "${sorted_mark}" | awk '/^SN/ {print \$2, \$3}' > "${sample_id}-markdup-samtools-stats.txt"
    
    ## define a Bash variable from the channel value
    sample_id_bash="${sample_id}"
    
    ## strip the literal suffix '_sorted'
    clean_id="\${sample_id_bash%_sorted}" 
    
    printf "SampleID\ttotal_depth\tpanel_depth\n%s\t%s\t%s\n" \
       "\${clean_id}" "\${total_depth}" "\${panel_depth}" \
            > "${sample_id}-depth.tsv"
    
    # ---------- read counts ---------- 
    total_reads=\$(samtools view -c -F 4,1024 -q 20    "${sorted_ch}")
    panel_reads=\$(samtools view -c -F 4,1024 -q 20 -L "${panel}" "${sorted_ch}")

    printf "SampleID\ttotal_reads\tpanel_reads\n%s\t%s\t%s\n" \
       "\${clean_id}" "\${total_reads}" "\${panel_reads}" \
          > "${sample_id}-read-counts.tsv"
    """
}
