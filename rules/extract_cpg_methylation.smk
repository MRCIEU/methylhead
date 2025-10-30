rule extract_cpg_methylation:
    input:
        bam=config['paths']['output'] + "/mark_duplicates/{sample}.markdup.bam",
        fasta=config['inputs']['fasta']
    output:
        config['paths']['output'] + "/extract_cpg_methylation/{sample}_CpG.methylKit.gz"
    params:
        bin_dir=config['containers']['wgbs']['bin_dir'],
        out_dir=config['paths']['output'] 
    resources:
        **get_resources("light")
    shell:
        apptainer_exec("wgbs",
        """
        {params.bin_dir}/MethylDackel extract --methylKit \
	    -o {params.out_dir}/extract_cpg_methylation/{wildcards.sample} \
	    {input.fasta} {input.bam}
        gzip {params.out_dir}/extract_cpg_methylation/{wildcards.sample}_CpG.methylKit
        """)

