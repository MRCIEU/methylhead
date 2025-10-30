rule estimate_cell_counts:
    input:
        config['paths']['output'] + "/methylation_matrix/methylation-matrix.csv.gz"
    output:
        config['paths']['output'] + "/estimate_cell_counts/counts.csv"
    params:
        bin_dir=config['containers']['meth']['bin_dir'],
        scripts_dir=config['paths']['scripts'],
        cell_reference=config['inputs']['cell_reference']
    resources:
        **get_resources("light")
    shell:
        apptainer_exec("meth", 
        """
        {params.bin_dir}/Rscript --vanilla \
            {params.scripts_dir}/estimate-cell-counts.r \
            {params.cell_reference} {input} {output}
        """)

