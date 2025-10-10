process samtools_stats {

    input:
        tuple val(sample_id), path(aligned_bam)
        tuple val(sample_id), path(marked_bam)
        tuple val(sample_id), path(sorted_bam)
        tuple val(sample_id), path(sorted_bam_bai)
        path panel

    output:
        tuple val(sample_id), path("${sample_id}-samtools-stats.txt")        , emit: aligned_bam_stats
        tuple val(sample_id), path("${sample_id}-markdup-samtools-stats.txt"), emit: marked_bam_stats
        tuple val(sample_id), path("${sample_id}-read-counts.tsv")           , emit: read_counts_file

    script:
    """
    # ---------- basic QC ----------
    samtools stats "${aligned_bam}" | awk '/^SN/ {print \$2, \$3}' > "${sample_id}-samtools-stats.txt"
    samtools stats "${marked_bam}" | awk '/^SN/ {print \$2, \$3}' > "${sample_id}-markdup-samtools-stats.txt"
    
    ## define a Bash variable from the channel value
    sample_id_bash="${sample_id}"
    
    ## strip the literal suffix '_sorted'
    clean_id="\${sample_id_bash%_sorted}" 
    
    # ---------- read counts ---------- 
    total_reads=\$(samtools view -c -F 4,1024 -q 20    "${sorted_bam}")
    panel_reads=\$(samtools view -c -F 4,1024 -q 20 -L "${panel}" "${sorted_bam}")

    printf "SampleID\ttotal_reads\tpanel_reads\n%s\t%s\t%s\n" \
       "\${clean_id}" "\${total_reads}" "\${panel_reads}" \
          > "${sample_id}-read-counts.tsv"
    """
}
