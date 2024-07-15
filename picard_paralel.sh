#!/bin/bash

nextflow run main.nf \
-c nextflow.config \
-resume \
--reads "" \
--pipeline picard \
--intervals "" \
--genome_folder "" 
