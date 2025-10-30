rule fastqc:
    input:
        r1=config['paths']['output'] + "/merge_lanes/{sample}_R1.fastq.gz",
        r2=config['paths']['output'] + "/merge_lanes/{sample}_R2.fastq.gz"
    output:
        html1=config['paths']['output'] + "/fastqc/{sample}_R1_fastqc.html",
        html2=config['paths']['output'] + "/fastqc/{sample}_R2_fastqc.html",
        zip1=config['paths']['output'] + "/fastqc/{sample}_R1_fastqc.zip",
        zip2=config['paths']['output'] + "/fastqc/{sample}_R2_fastqc.zip"
    params:
        bin_dir=config['containers']['wgbs']['bin_dir'],
        out_dir=config['paths']['output']
    resources:
        **get_resources("light")
    shell:
        apptainer_exec("wgbs",
        """
        {params.bin_dir}/fastqc -o {params.out_dir}/fastqc {input.r1} {input.r2}
        """)

