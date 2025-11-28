rule qc_report:
    input:
        meth=config['paths']['output'] + "/methylation/matrix.csv.gz",
        coverage=config['paths']['output'] + "/methylation/coverage.csv.gz",
        cell_counts=config['paths']['output'] + "/cell_counts/matrix.csv",
        filtered_reads=config['paths']['output'] + "/multiqc/multiqc_data/cutadapt_filtered_reads_plot.txt",
        raw_reads=config['paths']['output'] + "/samtools_stats/read-counts.tsv",
        hsmetrics=config['paths']['output'] + "/multiqc/multiqc_data/picard_hsmetrics_table.txt",
        panel=config['inputs']['panel']
    output:
        config['paths']['output'] + "/qc_report/qc-report.html"
    params:
        bin_dir=config['containers']['qc']['bin_dir'],
        scripts_dir=config['paths']['scripts'],
        out_dir=config['paths']['output']
    resources:
        **get_resources("light")
    shell:
        apptainer_exec("qc",
        """
	ln -sf {params.scripts_dir}/qc.qmd {params.out_dir}/qc_report/qc.qmd
        cd {params.out_dir}/qc_report
        {params.bin_dir}/quarto render qc.qmd \
            --execute-dir . \
            --to html \
            --output qc-report.html \
            --output-dir . \
            -P meth={input.meth} \
            -P coverage={input.coverage} \
            -P cell_counts={input.cell_counts} \
            -P filtered_reads={input.filtered_reads} \
            -P raw_reads={input.raw_reads} \
            -P hs_metrics={input.hsmetrics} \
            -P panel={input.panel}
        """)
