rule mbias:
    input:
        bam=config['paths']['output'] + "/mark_duplicates/{sample}.markdup.bam",
        fasta=config['inputs']['fasta']
    output:
        ot=config['paths']['output'] + "/mbias/{sample}_OT.svg",
        ob=config['paths']['output'] + "/mbias/{sample}_OB.svg"
    params:
        bin_dir=config['containers']['wgbs']['bin_dir'],
        out_dir=config['paths']['output']
    resources:
        **get_resources("light")
    shell:
        apptainer_exec("wgbs",
        """
        {params.bin_dir}/MethylDackel mbias \
	    {input.fasta} {input.bam} \
	    {params.out_dir}/mbias/{wildcards.sample} \
            --nOT 0,0,0,98 --nOB 0,0,3,0
        """)
