rule meth_bedgraph:
    input:
        config['paths']['output'] + "/methyldackel_bedgraph/{sample}_CpG.bedGraph.gz"
    output:
        config['paths']['output'] + "/meth_bedgraph/{sample}_cpg_meth.bedgraph.gz"
    params:
        scripts_dir=config['paths']['scripts']
    resources:
        **get_resources("light")
    shell:
        apptainer_exec("wgbs",
        """
        gunzip -c {input} | \
            {params.scripts_dir}/reformat-methyldackel-bedgraph.sh | \
            gzip -c > {output}
        """)

