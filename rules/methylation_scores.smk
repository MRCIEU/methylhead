rule methylation_scores:
    input:
        config['paths']['output'] + "/illumina/matrix.csv.gz"
    output:
        scores=config['paths']['output'] + "/methylation_scores/matrix.csv",
        sites=config['paths']['output'] + "/methylation_scores/sites.csv"
    params:
        bin_dir=config['containers']['meth']['bin_dir'],
        scripts_dir=config['paths']['scripts']
    resources:
        **get_resources("light")
    shell:
        apptainer_exec("meth",
        """
        {params.bin_dir}/Rscript --vanilla \
	    {params.scripts_dir}/methylation-scores.r \
	    {input} {output.scores} {output.sites}
        """)

