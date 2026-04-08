rule test_associations:
    input:
        phenotype=config['inputs']['phenotype'],
        models=config['inputs']['models'],
        meth=config['paths']['output'] + "/methylation/matrix.csv.gz",
        entropy=config['paths']['output'] + "/methyl_entropy/matrix.csv.gz",
        flips=config['paths']['output'] + "/methyl_entropy/flips.csv.gz",
        illumina=config['paths']['output'] + "/illumina/matrix.csv.gz",
        coverage=config['paths']['output'] + "/methylation/coverage.csv.gz",
        scores=config['paths']['output'] + "/methylation_scores/matrix.csv",
        methylseqlm=config['paths']['output'] + "/methylseqlm/matrix.csv.gz",
        cell_counts=config['paths']['output'] + "/cell_counts/matrix.csv",
        reads=config['paths']['output'] + "/multiqc/multiqc_data/cutadapt_filtered_reads_plot.txt"
    output:
        directory(config['paths']['output'] + "/test_associations/summary-stats/")
    params:
        bin_dir=config['containers']['meth']['bin_dir'],
        scripts_dir=config['paths']['scripts']
    resources:
        **get_resources("medium")
    shell:
        apptainer_exec("meth",
        """
        {params.bin_dir}/Rscript {params.scripts_dir}/test-associations.r \
            {input.phenotype} \
	    {input.models} \
	    {input.meth} \
	    {input.entropy} \
	    {input.flips} \
	    {input.illumina} \
	    {input.scores} \
	    {input.methylseqlm} \
	    {input.cell_counts} \
	    {input.coverage} \
	    {input.reads} \
	    {output}
        """)
