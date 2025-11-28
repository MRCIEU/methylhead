rule methylseqlm:
    input:
        phenotype=config['inputs']['phenotype'],
        model=config['inputs']['region_model'],
        meth=config['paths']['output'] + "/methylation/matrix.csv.gz",
        coverage=config['paths']['output'] + "/methylation/coverage.csv.gz",
        cell_counts=config['paths']['output'] + "/cell_counts/matrix.csv",
        reads=config['paths']['output'] + "/multiqc/multiqc_data/cutadapt_filtered_reads_plot.txt"
    output:
        config['paths']['output'] + "/methylseqlm/matrix.csv.gz"
    params:
        bin_dir=config['containers']['meth']['bin_dir'],
        scripts_dir=config['paths']['scripts']
    resources:
        **get_resources("medium")
    shell:
        apptainer_exec("meth",
        """
        {params.bin_dir}/Rscript {params.scripts_dir}/methylseqlm-matrix.r \
	    {params.scripts_dir} \
            {input.phenotype} \
	    {input.model} \
	    {input.meth} \
	    {input.cell_counts} \
	    {input.coverage} \
	    {input.reads} \
	    {output} \
	    {resources.cpus}
        """)


