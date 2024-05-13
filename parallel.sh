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


WORKFLOW="main.nf"
CONFIG=$2
DATA_DIR="/user/work/ag24712/Next_1/" 
REFERENCE_DIR="/user/work/ag24712/reference/hg38/"  
RESULTS_DIR="results" 

DATA_SUBDIR=data$SLURM_ARRAY_TASK_ID

start_time=$(date +"%s")

#conda init ## Required to make conda happy on the nodes
#source ~/.bashrc ## Required to load what conda init just did
#conda activate Bismark ## see for more details https://dsbristol.github.io/dst/coursebook/appendix5-bluecrystal.html

## run with -resume to avoid repeating previous analyses
nextflow run ${WORKFLOW} \
    --reads "${DATA_DIR}/${DATA_SUBDIR}/*_{1,2}.fastq.gz" \
    --t_param 8 \
    --memory_param 10000 \
    --multicore 8 \
    --cores 4 \
    --genome_folder ${REFERENCE_DIR} \
    --outdir ${RESULTS_DIR}

cp .nextflow.log ${RESULTS_DIR}/nextflow.log

cd /user/work/ag24712/Next_1/results/results/

bismark2report

bismark2summary

cd ..
cd ..

multiqc .


