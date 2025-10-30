rule illumina_matrix:
    input:
        config['paths']['output'] + "/methylation_matrix/methylation-matrix.csv.gz"
    output:
        config['paths']['output'] + "/illumina_matrix/illumina-matrix.csv.gz"
    params:
        bin_dir=config['containers']['meth']['bin_dir'],
        scripts_dir=config['paths']['scripts'],
        assembly=config['parameters']['assembly']
    resources:
        **get_resources("light")
    shell:
        apptainer_exec("meth",
        """
        {params.bin_dir}/Rscript --vanilla \
	    {params.scripts_dir}/illumina-matrix.r \
            {input} {output} {params.assembly}
        """)
