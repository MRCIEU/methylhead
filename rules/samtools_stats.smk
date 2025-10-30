rule samtools_stats:
    input:
        aligned=config['paths']['output'] + "/bwa_align/{sample}.bam",
        marked=config['paths']['output'] + "/mark_duplicates/{sample}.markdup.bam",
        sorted=config['paths']['output'] + "/sort_bam/{sample}_sorted.bam",
        bai=config['paths']['output'] + "/sort_bam/{sample}_sorted.bam.bai",
        panel=config['paths']['output'] + "/panel_intervals/panel.bed"
    output:
        aligned_stats=config['paths']['output'] + "/samtools_stats/{sample}-samtools-stats.txt",
        marked_stats=config['paths']['output'] + "/samtools_stats/{sample}-markdup-samtools-stats.txt",
        read_counts=config['paths']['output'] + "/samtools_stats/{sample}-read-counts.tsv"
    params:
        bin_dir=config['containers']['wgbs']['bin_dir'],
        scripts_dir=config['paths']['scripts']
    resources:
        **get_resources("light")
    shell:
        apptainer_exec("wgbs",
        """
        {params.bin_dir}/samtools stats {input.aligned} \
	    | {params.scripts_dir}/samtools_stats2csv.sh \
	    > {output.aligned_stats}
        {params.bin_dir}/samtools stats {input.marked} \
	    | {params.scripts_dir}/samtools_stats2csv.sh \
	    > {output.marked_stats}

        total_reads=\\$({params.bin_dir}/samtools view -c -F 4,1024 -q 20 {input.sorted})
        panel_reads=\\$({params.bin_dir}/samtools view -c -F 4,1024 -q 20 -L {input.panel} {input.sorted})
        
        printf 'SampleID\\ttotal_reads\\tpanel_reads\\n%s\\t%s\\t%s\\n' \
            {wildcards.sample} \\$total_reads \\$panel_reads \
            > {output.read_counts}
        """)
