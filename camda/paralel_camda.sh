#!/bin/bash

nextflow run main.nf \
-c nextflow.config \
-resume \
--reads "" \
--genome_folder "" \
--samtools_path ""
