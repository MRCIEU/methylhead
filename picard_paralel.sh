#!/bin/bash

nextflow run main.nf \
-c nextflow.config \
-resume \
--reads "/user/work/ag24712/fastq_files/*_R{1,2}*.fastq.gz" \
--pipeline picard \
--intervals "/user/work/ag24712/reference/hg19/Twist/covered_targets_Twist_Methylome_hg19_annotated_collapsed_uniq_final" --genome_folder "/user/work/ag24712/reference/hg19/Twist/bwa_index/hg19.fa" \
--cores 24
