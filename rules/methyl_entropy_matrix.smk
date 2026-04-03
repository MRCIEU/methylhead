rule collect_entropy_files:
    input:
        expand("{out_dir}/extract_methyl_entropy/{sample}.csv.gz",
            sample=all_samples,
            out_dir=config['paths']['output'])
    output:
        files=config['paths']['output'] + "/methyl_entropy_files.csv"
    run:
        with open(output.files, "w") as f:
            f.write("sample_id,filename\n")
            for fn in input:
                id=os.path.basename(fn).replace(".csv.gz","")
                f.write(f"{id},{fn}\n")

rule collect_entropy_regions:
    input:
        config['paths']['output'] + "/methyl_entropy_files.csv"
    output:
        config['paths']['output'] + "/methyl_entropy_regions.csv"
    resources:
        **get_resources("medium")
    shell:
        apptainer_exec("meth",
        """
        csvtk concat @{input} \
        | csvtk cut -t -f 1,2,3 \
        | csvtk freq -t -f 1,2,3 \
        > {output}
        """)

rule methyl_entropy_matrix:
    input:
        files=config['paths']['output'] + "/methyl_entropy_files.csv",
        regions=config['paths']['output'] + "/methyl_entropy_regions.csv"
    output:
        entropy=config['paths']['output'] + "/methyl_entropy/matrix.csv.gz",
        meth   =config['paths']['output'] + "/methyl_entropy/meth.csv.gz",
        sites  =config['paths']['output'] + "/methyl_entropy/sites.csv.gz"
    params:
        bin_dir=config['containers']['meth']['bin_dir'],
        scripts_dir=config['paths']['scripts']
    resources:
        **get_resources("medium")
    shell:
        apptainer_exec("meth",
        """
        {params.bin_dir}/Rscript \
            {params.scripts_dir}/methyl-entropy-matrix.r \
            {input.files} {input.regions}  \
            {output.entropy} {output.meth} {output.sites} \
            {resources.cpus}
        """)
