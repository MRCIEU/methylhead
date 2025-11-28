rule methyldackel_bedgraph:
    input:
        bam=config['paths']['output'] + "/mark_duplicates/{sample}.markdup.bam",
        fasta=config['inputs']['fasta']
    output:
        temp(config['paths']['output'] + "/methyldackel_bedgraph/{sample}_CpG.bedGraph.gz")
    params:
        bin_dir=config['containers']['wgbs']['bin_dir'],
        out_dir=config['paths']['output']
    resources:
        **get_resources("light")
    shell:
        apptainer_exec("wgbs", 
        """
        {params.bin_dir}/MethylDackel extract \
	    -o {params.out_dir}/methyldackel_bedgraph/{wildcards.sample} \
            --minDepth 10 --maxVariantFrac 0.25 \
            --nOT 0,0,0,98 --nOB 0,0,3,0 --mergeContext \
            {input.fasta} {input.bam} --keepDupes
        gzip {params.out_dir}/methyldackel_bedgraph/{wildcards.sample}_CpG.bedGraph
        """)
