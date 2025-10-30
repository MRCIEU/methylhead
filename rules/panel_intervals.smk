rule panel_intervals:
    input:
        panel=config['inputs']['panel'],
        fasta=config['inputs']['fasta'],
        dict=config['inputs']['fasta'] + ".dict"
    output:
        interval=config['paths']['output'] + "/panel_intervals/interval_file",
        bed=config['paths']['output'] + "/panel_intervals/panel.bed"
    params:
        bin_dir=config['containers']['wgbs']['bin_dir'],
        scripts_dir=config['paths']['scripts'],
        out_dir=config['paths']['output']
    resources:
        **get_resources("light")
    shell:
        apptainer_exec("wgbs",
        """
	bash {params.scripts_dir}/csv2bed.sh \
	     {input.panel} \
	     {output.bed}
	{params.bin_dir}/picard BedToIntervalList \
    	    I={output.bed} \
    	    O={output.interval} \
    	    SD={input.fasta}
	""")
