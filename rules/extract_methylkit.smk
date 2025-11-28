rule extract_methylkit:
    input:
        bam=config['paths']['output'] + "/mark_duplicates/{sample}.markdup.bam",
        fasta=config['inputs']['fasta']
    output:
        config['paths']['output'] + "/extract_methylkit/{sample}_CpG.methylKit.gz"
    params:
        bin_dir=config['containers']['wgbs']['bin_dir'],
        out_dir=config['paths']['output'] 
    resources:
        **get_resources("light")
    shell:
        apptainer_exec("wgbs",
        """
        {params.bin_dir}/MethylDackel extract --methylKit \
	    -o {params.out_dir}/extract_methylkit/{wildcards.sample} \
	    {input.fasta} {input.bam}
        gzip {params.out_dir}/extract_methylkit/{wildcards.sample}_CpG.methylKit
        """)

