rule mark_duplicates:
    input:
        config['paths']['output'] + "/sort_bam/{sample}_sorted.bam"
    output:
        bam=config['paths']['output'] + "/mark_duplicates/{sample}.markdup.bam",
        bai=config['paths']['output'] + "/mark_duplicates/{sample}.markdup.bam.bai",
        metrics=config['paths']['output'] + "/mark_duplicates/{sample}.markdup_metrics.txt"
    params:
        bin_dir=config['containers']['wgbs']['bin_dir']
    resources:
        **get_resources("medium")
    shell:
        apptainer_exec("wgbs",
        """
        {params.bin_dir}/picard MarkDuplicates \
            INPUT={input} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metrics} \
            REMOVE_DUPLICATES=false \
            ASSUME_SORT_ORDER=coordinate \
            OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500        
        {params.bin_dir}/samtools index -@ {resources.cpus} {output.bam}
        """)
