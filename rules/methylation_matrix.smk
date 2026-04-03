rule collect_meth_files:
    input:
        expand("{out_dir}/extract_methylkit/{sample}_CpG.methylKit.gz",
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

rule collect_methyl_sites:
    input:
        config['paths']['output'] + "/meth_files.csv" 
    output:
        config['paths']['output'] + "/meth_sites.csv"
    resources:
        **get_resources("medium")
    shell:
        apptainer_exec("meth",
        """
        csvtk concat @{input} \
        | csvtk cut -t -f 2,3 \
        | csvtk freq -t -f 1,2 \
        > {output}
        """)

rule methylation_matrix:
    input:
        files=config['paths']['output'] + "/meth_files.csv",
        regions=config['paths']['output'] + "/meth_sites.csv"
    output:
        meth=config['paths']['output'] + "/methylation/matrix.csv.gz",
        coverage=config['paths']['output'] + "/methylation/coverage.csv.gz"
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
            {input.files} {input.regions} \
            {output.meth} {output.coverage} \
            {resources.cpus}
        """)

