#!/bin/bash

nextflow run main.nf \
-c nextflow.config \
-resume \
--reads "" \
--pipeline bismark \
--genome_folder "" 
