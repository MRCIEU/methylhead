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


nextflow run ${WORKFLOW} \
    --reads "${DATA_DIR}/${DATA_SUBDIR}/*_{1,2}.fastq.gz" \
    --t_param 8 \
    --memory_param 10000 \
    --multicore 8 \
    --cores 4 \
    --genome_folder ${REFERENCE_DIR} \
    --outdir ${RESULTS_DIR}


end_time=$(date +"%s")


cpu_time=$(sacct -j $SLURM_JOBID --format=JobID,TotalCPU | grep $SLURM_JOBID | awk '{print $2}')


total_time=$(($end_time-$start_time))


echo "CPU Zamani: $cpu_time"
echo "Toplam Zaman: $total_time saniye"


cp .nextflow.log ${RESULTS_DIR}/nextflow.log


echo "CPU Zamani: $cpu_time" >> ${RESULTS_DIR}/execution.log
echo "Toplam Zaman: $total_time saniye" >> ${RESULTS_DIR}/execution.log

cd /user/work/ag24712/Next_1/results/results/

bismark2report

bismark2summary

cd ..
cd ..

multiqc .


