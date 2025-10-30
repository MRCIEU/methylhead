## assemble lists of fastq file paths for each sample
fastq_dir = dirname(config["inputs"]["samplesheet"])
samples = pandas.read_csv(config["inputs"]["samplesheet"])
samples["read1"] = [normalize_path(path, fastq_dir) for path in samples["read1"]]
samples["read2"] = [normalize_path(path, fastq_dir) for path in samples["read2"]]
read1_paths = samples.groupby("sample_id")["read1"].apply(list).to_dict()
read2_paths = samples.groupby("sample_id")["read2"].apply(list).to_dict()

rule merge_lanes:
    input:
        r1=lambda wildcards: read1_paths[wildcards.sample],
        r2=lambda wildcards: read2_paths[wildcards.sample]
    output:
        r1=temp(config['paths']['output'] + "/merge_lanes/{sample}_R1.fastq.gz"),
        r2=temp(config['paths']['output'] + "/merge_lanes/{sample}_R2.fastq.gz")
    resources:
        **get_resources("light")
    shell:
        apptainer_exec("wgbs", 
	"""
        cat {input.r1} > {output.r1}
        cat {input.r2} > {output.r2}
        """)
