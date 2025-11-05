#!/bin/bash

CONFIG=config.yml
MAX_CORES=16

snakemake --snakefile ../Snakefile \
  --configfile $CONFIG --cores $MAX_CORES \
  --printshellcmds --rerun-incomplete \
  --rerun-triggers mtime \
  --show-failed-logs --verbose 

