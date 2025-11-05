rule collect_camda_files:
    input:
        files=expand("{out_dir}/camda/{sample}_CpG_CAMDA.tsv.gz",
	    sample=all_samples,
	    out_dir=config['paths']['output'])
    output:
        files=config['paths']['output'] + "/camda_files.csv"
    run:
        with open(output.files, "w") as f:
            f.write("sample_id,filename\n")
            for fn in input:
                id=os.path.basename(fn).replace("_CpG_CAMDA.tsv.gz","")
                f.write(f"{id},{fn}\n")

rule camda_matrix:
    input:
        config['paths']['output'] + "/camda_files.csv"
    output:
        config['paths']['output'] + "/camda_matrix/camda-matrix.csv.gz"
    params:
        bin_dir=config['containers']['meth']['bin_dir'],
        scripts_dir=config['paths']['scripts']
    resources:
        **get_resources("medium")
    shell:
        apptainer_exec("meth",
        """
        {params.bin_dir}/Rscript --vanilla {params.scripts_dir}/camda-matrix.r \
	    {input}  {params.scripts_dir} {output} 
        """)

