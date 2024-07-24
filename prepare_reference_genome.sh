#!/bin/bash

# Check if the user provided the GENOMES path
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <path_to_converted_genome_index>"
    exit 1
fi

# Define the genome index path from user input
GENOMES=$1

# Singularity image path
SINGULARITY_IMAGE="dnam_cancer_pipeline_latest.sif"

# Create the GENOMES directory if it does not exist
mkdir -p ${GENOMES}

# Download the hg19 reference genome if it does not exist
if [ ! -f ${GENOMES}/hg19.fa ]; then
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz -P ${GENOMES}
    gunzip ${GENOMES}/hg19.fa.gz 
fi

# Index the reference genome using Singularity if index does not exist
if [ ! -f ${GENOMES}/hg19.fa.fai ]; then
    singularity exec ${SINGULARITY_IMAGE} samtools faidx ${GENOMES}/hg19.fa
fi

# Prepare the genome with bismark using Singularity if Bismark index does not exist
if [ ! -d ${GENOMES}/Bisulfite_Genome ]; then
    singularity exec ${SINGULARITY_IMAGE} bismark_genome_preparation --path_to_aligner /opt/bismark/bismark_genome_preparation --verbose ${GENOMES}
fi

# Download the BED file if it does not exist
if [ ! -f ${GENOMES}/covered_targets_Twist_Methylome_hg19_annotated_collapsed_final.bed ]; then
    wget https://www.twistbioscience.com/sites/default/files/resources/2022-06/covered_targets_Twist_Methylome_hg19_annotated_collapsed_final.bed.zip -P ${GENOMES}
    unzip ${GENOMES}/covered_targets_Twist_Methylome_hg19_annotated_collapsed_final.bed.zip -d ${GENOMES}
fi

# Create the sequence dictionary using Singularity if it does not exist
if [ ! -f ${GENOMES}/hg19.fa.dict ]; then
    singularity exec ${SINGULARITY_IMAGE} picard CreateSequenceDictionary REFERENCE=${GENOMES}/hg19.fa OUTPUT=${GENOMES}/hg19.fa.dict
fi

# Convert BED to interval list using Singularity if it does not exist
if [ ! -f ${GENOMES}/covered_targets_Twist_Methylome_hg19_annotated_collapsed_final.interval_list ]; then
    singularity exec ${SINGULARITY_IMAGE} picard BedToIntervalList \
        I=${GENOMES}/covered_targets_Twist_Methylome_hg19_annotated_collapsed_final.bed \
        O=${GENOMES}/covered_targets_Twist_Methylome_hg19_annotated_collapsed_final.interval_list \
        SD=${GENOMES}/hg19.fa.dict
fi

# Index the reference genome using BWA-MEM with bwameth.py if BWA-MEM index does not exist
if [ ! -f ${GENOMES}/hg19.fa.bwt ]; then
    singularity exec ${SINGULARITY_IMAGE} bwameth.py index ${GENOMES}/hg19.fa
fi

echo "Reference genome preparation completed."
