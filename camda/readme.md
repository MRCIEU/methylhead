Parameters


--reads: Specifies the path to the input FASTQ files.


--genome_folder: Specifies the path to the reference genome file.


--samtools_path: Specifies the path to the samtools.



** This scripts downloaded from https://github.com/JiejunShi/CAMDA/tree/master



Usage
```
 nextflow main.nf --reads "~/path/*_R{1,2}*.fastq.gz" -resume --outdir CAMDA --genome_folder "~/path/hg19.fa" --samtools_path "~/path/"
```
