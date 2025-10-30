rule bwa_align:
    input:
        r1=config['paths']['output'] + "/trim/{sample}_R1_val_1.fq.gz",
        r2=config['paths']['output'] + "/trim/{sample}_R2_val_2.fq.gz",
        fasta=config['inputs']['fasta'],
        bwaindex=config['inputs']['fasta'] + ".bwameth.c2t"
    output:
        temp(config['paths']['output'] + "/bwa_align/{sample}.bam")
    params:
        bin_dir=config['containers']['wgbs']['bin_dir']
    resources:
        **get_resources("heavy")
    shell:
        apptainer_exec("wgbs",
        """
	{params.bin_dir}/bwameth.py \
	    --reference {input.fasta} {input.r1} {input.r2} -t {threads} | \
	    {params.bin_dir}/samtools view -b - > {output}
        """)
