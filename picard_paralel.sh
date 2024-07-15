#!/bin/bash

nextflow run main.nf \
-c nextflow.config \
-resume \
--data "" \
--pipeline picard \
--intervals "" \
--genome_folder "" 
