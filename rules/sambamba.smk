rule sambamba:
    input:
        bam=config['paths']['output'] + "/bwa_align/{sample}.bam",
        fasta=config['inputs']['fasta']
    output:
        temp(config['paths']['output'] + "/sambamba/{sample}.bam")
    params:
        bin_dir=config['containers']['wgbs']['bin_dir']
    resources:
        **get_resources("medium")
    shell:
        apptainer_exec("wgbs",
        """
        {params.bin_dir}/sambamba view -h \
	    -t {resources.cpus} -T {input.fasta} \
            --filter='not secondary_alignment and not failed_quality_control and not supplementary and proper_pair and mapping_quality > 0' \
            -f bam -l 0 {input.bam} -o {output}
        """)
