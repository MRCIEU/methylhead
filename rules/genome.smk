rule genome:
    output:
        fasta=config['inputs']['fasta'],
        samindex=config['inputs']['fasta'] + ".fai",
        dictionary=config['inputs']['fasta'] + ".dict",
        bwaindex=config['inputs']['fasta'] + ".bwameth.c2t"
    params:
        scripts_dir=config['paths']['scripts'],
        assembly=config['parameters']['assembly'],
        genome_dir=config['paths']['genome']
    resources:
        **get_resources("heavy")
    shell:
        apptainer_exec("wgbs",
        """
        bash {params.scripts_dir}/create-reference.sh \
            {params.assembly} {params.genome_dir}
        """)
