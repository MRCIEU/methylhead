rule collect_hs_metrics:
    input:
        bam=config['paths']['output'] + "/mark_duplicates/{sample}.markdup.bam",
        fasta=config['inputs']['fasta'],
        intervals=config['paths']['output'] + "/panel_intervals/interval_file"
    output:
        metrics=config['paths']['output'] + "/collect_hs_metrics/{sample}_coverage_metrics.txt",
        coverage=config['paths']['output'] + "/collect_hs_metrics/{sample}_coverage"
    params:
        bin_dir=config['containers']['wgbs']['bin_dir']
    resources:
        **get_resources("medium")
    shell:
        apptainer_exec("wgbs", 
        """
        {params.bin_dir}/picard CollectHsMetrics \
            I={input.bam} \
            O={output.metrics} \
	    R={input.fasta}  \
            BAIT_INTERVALS={input.intervals} \
            TARGET_INTERVALS={input.intervals} \
            MINIMUM_MAPPING_QUALITY=20 \
            COVERAGE_CAP=200 \
            PER_TARGET_COVERAGE={output.coverage} \
            NEAR_DISTANCE=250
        """)
