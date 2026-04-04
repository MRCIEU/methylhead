rule extract_methyl_entropy:
    input:
        bam=config['paths']['output'] + "/mark_duplicates/{sample}.markdup.bam",
        fasta=config['inputs']['fasta']
    output:
        config['paths']['output'] + "/extract_methyl_entropy/{sample}.csv.gz"
    params:
        scripts_dir=config['paths']['scripts'],
        region_size=config['parameters']['region_size']
    resources:
        **get_resources("medium")
    shell:
        apptainer_exec("meth",
        """
        /opt/venv/bin/python {params.scripts_dir}/extract-clonal-flips.py \
             --size {params.region_size} {input.bam} {input.fasta} {output} 
         """)

