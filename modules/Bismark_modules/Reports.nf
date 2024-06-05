process Reports {

        input:
        
        val alignment_t
        val deduplication_t
        val methylation_t
        val methylation_t2
        val bam2nuc

        
        publishDir "${params.outdir}/Reports" , mode: 'copy'
      
        output:
        
        path("*.html"), emit: reports
        
        script:
        """
       bismark2report --alignment_report "${alignment_t}" --dedup_report "${deduplication_t}" --splitting_report "${methylation_t}" --mbias_report "${methylation_t2}" --nucleotide_report "${bam2nuc}"      
        """
}


