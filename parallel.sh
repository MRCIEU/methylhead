#!/bin/bash

#SBATCH --job-name=nextflow_parallel
#SBATCH --account=SSCM009461
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#SBATCH --time=0:60:00
#SBATCH --mem=100G
#SBATCH --array=1-3

conda init
conda activate Bismark

source config.sh

DATA_SUBDIR=data$SLURM_ARRAY_TASK_ID

nextflow run ${WORKFLOW} \
    --reads "${DATA_DIR}/${DATA_SUBDIR}/*_{1,2}.fastq.gz" \
    --t_param 8 \
    --memory_param 10000 \
    --multicore 8 \
    --cores 4 \
    --genome_folder ${REFERENCE_DIR} \
    --outdir ${RESULTS_DIR}
