rule bsmap_align:
    input:
        r1=config['paths']['output'] + "/trim/{sample}_R1_val_1.fq.gz",
        r2=config['paths']['output'] + "/trim/{sample}_R2_val_2.fq.gz",
        fasta=config['inputs']['fasta']
    output:
        temp(config['paths']['output'] + "/bsmap_align/{sample}.bam")
    params:
        bin_dir=config['containers']['wgbs']['bin_dir']
    resources:
        **get_resources("heavy")
    shell:
        apptainer_exec("wgbs",
        """
        {params.bin_dir}/bsmap \
	    -d {input.fasta} \
	    -a {input.r1} -b {input.r2} \
            -R -p {resources.cpus} -o {output}
        """)
