rule sort_bam:
    input:
        config['paths']['output'] + "/sambamba/{sample}.bam"
    output:
        bam=temp(config['paths']['output'] + "/sort_bam/{sample}_sorted.bam"),
        bai=temp(config['paths']['output'] + "/sort_bam/{sample}_sorted.bam.bai")
    params:
        bin_dir=config['containers']['wgbs']['bin_dir']
    resources:
        **get_resources("medium")
    shell:
        apptainer_exec("wgbs",
        """
	mem_per_cpu=\\$(({resources.mem_gb} / {resources.cpus}))G
        {params.bin_dir}/samtools sort \
	    -@ {resources.cpus} \
	    -m \\$mem_per_cpu \
	    -l 0 \
	    -o {output.bam} \
	    {input}
        {params.bin_dir}/samtools index \
	    -@ {resources.cpus} \
	    {output.bam} > {output.bai}
        """)
