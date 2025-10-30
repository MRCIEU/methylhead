rule dna_methylation_scores:
    input:
        config['paths']['output'] + "/illumina_matrix/illumina-matrix.csv.gz"
    output:
        scores=config['paths']['output'] + "/dna_methylation_scores/scores.csv",
        sites=config['paths']['output'] + "/dna_methylation_scores/sites.csv"
    params:
        bin_dir=config['containers']['meth']['bin_dir'],
        scripts_dir=config['paths']['scripts']
    resources:
        **get_resources("light")
    shell:
        apptainer_exec("meth",
        """
        {params.bin_dir}/Rscript --vanilla \
	    {params.scripts_dir}/dna-methylation-scores.r \
	    {input} {output.scores} {output.sites}
        """)

