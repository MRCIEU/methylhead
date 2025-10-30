rule collect_mm_metrics:
    input:
        bam=config['paths']['output'] + "/mark_duplicates/{sample}.markdup.bam",
        fasta=config['inputs']['fasta']
    output:
        alignment=config['paths']['output'] + "/collect_mm_metrics/{sample}.alignment_summary_metrics",
        gc_detail=config['paths']['output'] + "/collect_mm_metrics/{sample}.gc_bias.detail_metrics",
        gc_pdf=config['paths']['output'] + "/collect_mm_metrics/{sample}.gc_bias.pdf",
        gc_summary=config['paths']['output'] + "/collect_mm_metrics/{sample}.gc_bias.summary_metrics",
        insert_hist=config['paths']['output'] + "/collect_mm_metrics/{sample}.insert_size_histogram.pdf",
        insert_metrics=config['paths']['output'] + "/collect_mm_metrics/{sample}.insert_size_metrics"
    params:
        out_dir=config['paths']['output'],
        bin_dir=config['containers']['wgbs']['bin_dir']
    resources:
        **get_resources("medium")
    shell:
        apptainer_exec("wgbs", 
        """
        {params.bin_dir}/picard CollectMultipleMetrics \
            I={input.bam} \
            O={params.out_dir}/collect_mm_metrics/{wildcards.sample} \
	    R={input.fasta} \
            PROGRAM=null \
            PROGRAM=CollectGcBiasMetrics \
            PROGRAM=CollectInsertSizeMetrics \
            PROGRAM=CollectAlignmentSummaryMetrics
        """)
