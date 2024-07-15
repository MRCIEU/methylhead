#!/bin/bash

nextflow run main.nf \
-c nextflow.config \
-resume \
--data "" \
--pipeline bismark \
--genome_folder "" 
