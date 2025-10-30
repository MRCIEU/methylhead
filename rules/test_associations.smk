rule test_associations:
    input:
        phenotype=config['inputs']['phenotype'],
        models=config['inputs']['models'],
        meth=config['paths']['output'] + "/methylation_matrix/methylation-matrix.csv.gz",
        illumina=config['paths']['output'] + "/illumina_matrix/illumina-matrix.csv.gz",
        coverage=config['paths']['output'] + "/methylation_matrix/coverage-matrix.csv.gz",
        scores=config['paths']['output'] + "/dna_methylation_scores/scores.csv",
        camda=config['paths']['output'] + "/camda_matrix/camda-matrix.csv.gz",
        cell_counts=config['paths']['output'] + "/estimate_cell_counts/counts.csv",
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
	    {input.illumina} \
	    {input.scores} \
	    {input.camda} \
	    {input.cell_counts} \
	    {input.coverage} \
	    {input.reads} \
	    {output}
        """)
