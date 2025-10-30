#!/bin/bash

snakemake --snakefile ../Snakefile \
  --configfile config.yml --cores 16 \
  --printshellcmds --rerun-incomplete \
  --rerun-triggers mtime \
  --show-failed-logs --verbose 
