rule camda:
    input:
        bam=config['paths']['output'] + "/bsmap_align/{sample}.bam",
        fasta=config['inputs']['fasta']
    output:
        camda=config['paths']['output'] + "/camda/{sample}_CpG_CAMDA.tsv.gz",
        meth=config['paths']['output'] + "/camda/{sample}_CpG_MethRatio.tsv.gz"
    params:
        bin_dir=config['containers']['wgbs']['bin_dir'],
        out_dir=config['paths']['output'],
        scripts_dir=config['paths']['scripts']
    resources:
        **get_resources("heavy")
    shell:
        apptainer_exec("wgbs", 
        """
        cd {params.out_dir}/camda
        {params.bin_dir}/python {params.scripts_dir}/camda.py \
	    CAMDA {input.bam} {input.fasta}	\
            -o {wildcards.sample} \
            -w {wildcards.sample} \
            -s {params.bin_dir}/samtools  \
            -X CG
        gzip {wildcards.sample}_CpG_CAMDA.tsv
        gzip {wildcards.sample}_CpG_MethRatio.tsv
        gzip {wildcards.sample}_CpG_CAMDA.wig
        gzip {wildcards.sample}_CpG_MethRatio.wig
        """)
