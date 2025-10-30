rule collect_meth_files:
    input:
        expand("{out_dir}/extract_cpg_methylation/{sample}_CpG.methylKit.gz",
	    sample=all_samples,
	    out_dir=config['paths']['output'])
    output:
        files=config['paths']['output'] + "/meth_files.csv"
    run:
        with open(output.files, "w") as f:
            f.write("sample_id,filename\n")
            for fn in input:
                id=os.path.basename(fn).replace("_CpG.methylKit.gz","")	    
                f.write(f"{id},{fn}\n")

rule methylation_matrix:
    input:
        config['paths']['output'] + "/meth_files.csv" 
    output:
        meth=config['paths']['output'] + "/methylation_matrix/methylation-matrix.csv.gz",
        coverage=config['paths']['output'] + "/methylation_matrix/coverage-matrix.csv.gz"
    params:
        bin_dir=config['containers']['meth']['bin_dir'],
        scripts_dir=config['paths']['scripts']
    resources:
        **get_resources("medium")
    shell:
        apptainer_exec("meth",
        """
	{params.bin_dir}/Rscript \
	    {params.scripts_dir}/methylation-matrix.r \
	    {input} {output.meth} {output.coverage} 
        """)
