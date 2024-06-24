#!/bin/bash

nextflow run main.nf \
-c nextflow.config \
-resume \
--reads "/user/work/ag24712/fastq_files/*_R{1,2}*.fastq.gz" \
--pipeline bismark \
--genome_folder "/user/work/ag24712/reference/hg19/Bismark/" \
--cores 24
