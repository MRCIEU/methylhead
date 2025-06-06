#!/usr/bin/env bash

set -euo pipefail

# Directory to store downloaded FASTQ files

TARGET_DIR="fastq-files"
mkdir -p "$TARGET_DIR"
cd "$TARGET_DIR"

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/004/SRR14580504/SRR14580504_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/004/SRR14580504/SRR14580504_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/005/SRR14580505/SRR14580505_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/005/SRR14580505/SRR14580505_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/008/SRR14580508/SRR14580508_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/008/SRR14580508/SRR14580508_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/009/SRR14580509/SRR14580509_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/009/SRR14580509/SRR14580509_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/011/SRR14580511/SRR14580511_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/011/SRR14580511/SRR14580511_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/014/SRR14580514/SRR14580514_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/014/SRR14580514/SRR14580514_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/016/SRR14580516/SRR14580516_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/016/SRR14580516/SRR14580516_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/021/SRR14580521/SRR14580521_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/021/SRR14580521/SRR14580521_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/022/SRR14580522/SRR14580522_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/022/SRR14580522/SRR14580522_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/030/SRR14580530/SRR14580530_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/030/SRR14580530/SRR14580530_2.fastq.gz

