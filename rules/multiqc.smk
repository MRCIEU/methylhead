rule multiqc:
    input:
        expand("{out_dir}/fastqc/{sample}_R{read}_fastqc.zip",
	    sample=all_samples,
	    read=[1,2],
	    out_dir=config['paths']['output']),
        expand("{out_dir}/trim/{sample}_R{read}.fastq.gz_trimming_report.txt",
	    sample=all_samples,
	    read=[1,2],
	    out_dir=config['paths']['output']),
        expand("{out_dir}/mark_duplicates/{sample}.markdup_metrics.txt",
	    sample=all_samples,
	    out_dir=config['paths']['output']),
        expand("{out_dir}/collect_hs_metrics/{sample}_coverage_metrics.txt",
	    sample=all_samples,
	    out_dir=config['paths']['output']),
        expand("{out_dir}/collect_mm_metrics/{sample}.alignment_summary_metrics",
	    sample=all_samples,
	    out_dir=config['paths']['output']),
        expand("{out_dir}/samtools_stats/{sample}-samtools-stats.txt",
	    sample=all_samples,
	    out_dir=config['paths']['output'])
    output:
        config['paths']['output'] + "/multiqc/multiqc_report.html",
        config['paths']['output'] + "/multiqc/multiqc_data/cutadapt_filtered_reads_plot.txt",
        config['paths']['output'] + "/multiqc/multiqc_data/picard_hsmetrics_table.txt"
    params:
        bin_dir=config['containers']['wgbs']['bin_dir'],
        out_dir=config['paths']['output']
    resources:
        **get_resources("light")
    shell:
        apptainer_exec("wgbs",
        """
        cd {params.out_dir}
        {params.bin_dir}/multiqc . -o multiqc/ --force
        """)
