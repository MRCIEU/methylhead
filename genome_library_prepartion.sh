#!/bin/bash

# Shell script to prepare reference genome

# Define the genome index path
GENOMES=[path_to_converted_genome_index]

# Define the path to the aligner
ALIGNER_PATH="~/envs/Bismark/bin"

# Download the hg19 reference genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

# Move the downloaded genome to the GENOMES directory
mv hg19.fa.gz ${GENOMES}

# Prepare the genome with bismark
bismark_genome_preparation --path_to_aligner ${ALIGNER_PATH} --verbose ${GENOMES}

# Download the BED file
wget https://www.twistbioscience.com/sites/default/files/resources/2022-06/covered_targets_Twist_Methylome_hg19_annotated_collapsed_final.bed.zip

# Move the downloaded BED file to the GENOMES directory
mv covered_targets_Twist_Methylome_hg19_annotated_collapsed_final.bed.zip ${GENOMES}

# Unzip the BED file
unzip ${GENOMES}/covered_targets_Twist_Methylome_hg19_annotated_collapsed_final.bed.zip -d ${GENOMES}

# Index the reference genome
samtools faidx ${GENOMES}/hg19.fa

# Create the sequence dictionary
picard CreateSequenceDictionary REFERENCE=${GENOMES}/hg19.fa OUTPUT=${GENOMES}/hg19.fa.dict

# Convert BED to interval list
picard BedToIntervalList \
    I=${GENOMES}/covered_targets_Twist_Methylome_hg19_annotated_collapsed_final.bed \
    O=${GENOMES}/covered_targets_Twist_Methylome_hg19_annotated_collapsed_final \
    SD=${GENOMES}/hg19.fa.dict

echo "Reference genome preparation completed. Running time is approximately 2 hours."
