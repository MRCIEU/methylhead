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

# Download the hg19 reference genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz -P ${GENOMES}

# Prepare the genome with bismark using Singularity
singularity exec ${SINGULARITY_IMAGE} bismark_genome_preparation --path_to_aligner /opt/conda/envs/bismark/bin --verbose ${GENOMES}

# Download the BED file
wget https://www.twistbioscience.com/sites/default/files/resources/2022-06/covered_targets_Twist_Methylome_hg19_annotated_collapsed_final.bed.zip -P ${GENOMES}

# Unzip the BED file
unzip ${GENOMES}/covered_targets_Twist_Methylome_hg19_annotated_collapsed_final.bed.zip -d ${GENOMES}

# Index the reference genome using Singularity
singularity exec ${SINGULARITY_IMAGE} samtools faidx ${GENOMES}/hg19.fa

# Create the sequence dictionary using Singularity
singularity exec ${SINGULARITY_IMAGE} picard CreateSequenceDictionary REFERENCE=${GENOMES}/hg19.fa OUTPUT=${GENOMES}/hg19.fa.dict

# Convert BED to interval list using Singularity
singularity exec ${SINGULARITY_IMAGE} picard BedToIntervalList \
    I=${GENOMES}/covered_targets_Twist_Methylome_hg19_annotated_collapsed_final.bed \
    O=${GENOMES}/covered_targets_Twist_Methylome_hg19_annotated_collapsed_final \
    SD=${GENOMES}/hg19.fa.dict

echo "Reference genome preparation completed. Running time is approximately 2 hours."
