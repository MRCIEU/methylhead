rule raw_read_counts:
    input:
        expand("{out_dir}/samtools_stats/{sample}-read-counts.tsv",
	    sample=all_samples,
	    out_dir=config['paths']['output'])
    output:
        file=config['paths']['output'] + "/samtools_stats/read-counts.tsv"
    resources:
        **get_resources("light")
    run:
        first = True
        with open(output.file, "w") as output:
            for fn in input:
                with open(fn, "r") as input:
                    lines = input.readlines()
                    if first:
                        output.writelines(lines)
                    else:
                        output.writelines(lines[1])
                first = False
