rule trim:
    input:
        r1=config['paths']['output'] + "/merge_lanes/{sample}_R1.fastq.gz",
        r2=config['paths']['output'] + "/merge_lanes/{sample}_R2.fastq.gz"
    output:
        r1=temp(config['paths']['output'] + "/trim/{sample}_R1_val_1.fq.gz"),
        r2=temp(config['paths']['output'] + "/trim/{sample}_R2_val_2.fq.gz"),
        report1=config['paths']['output'] + "/trim/{sample}_R1.fastq.gz_trimming_report.txt",
        report2=config['paths']['output'] + "/trim/{sample}_R2.fastq.gz_trimming_report.txt"
    params:
        bin_dir=config['containers']['wgbs']['bin_dir'],
        out_dir=config['paths']['output']
    resources:
        **get_resources("medium")
    shell:
        apptainer_exec("wgbs",
        """
	{params.bin_dir}/trim_galore \
	    --paired {input.r1} {input.r2} \
	    --gzip -o {params.out_dir}/trim
        """)
