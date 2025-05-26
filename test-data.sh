#!/usr/bin/env bash

set -euo pipefail

# Directory to store downloaded FASTQ files

TARGET_DIR="test-data"
mkdir -p "$TARGET_DIR"
cd "$TARGET_DIR"

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/003/SRR14580503/SRR14580503_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/003/SRR14580503/SRR14580503_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/004/SRR14580504/SRR14580504_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/004/SRR14580504/SRR14580504_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/005/SRR14580505/SRR14580505_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/005/SRR14580505/SRR14580505_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/006/SRR14580506/SRR14580506_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/006/SRR14580506/SRR14580506_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/008/SRR14580508/SRR14580508_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/008/SRR14580508/SRR14580508_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/009/SRR14580509/SRR14580509_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/009/SRR14580509/SRR14580509_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/010/SRR14580510/SRR14580510_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/010/SRR14580510/SRR14580510_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/011/SRR14580511/SRR14580511_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/011/SRR14580511/SRR14580511_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/014/SRR14580514/SRR14580514_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/014/SRR14580514/SRR14580514_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/016/SRR14580516/SRR14580516_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/016/SRR14580516/SRR14580516_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/020/SRR14580520/SRR14580520_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/020/SRR14580520/SRR14580520_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/021/SRR14580521/SRR14580521_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/021/SRR14580521/SRR14580521_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/022/SRR14580522/SRR14580522_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/022/SRR14580522/SRR14580522_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/024/SRR14580524/SRR14580524_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/024/SRR14580524/SRR14580524_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/026/SRR14580526/SRR14580526_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/026/SRR14580526/SRR14580526_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/027/SRR14580527/SRR14580527_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/027/SRR14580527/SRR14580527_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/028/SRR14580528/SRR14580528_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/028/SRR14580528/SRR14580528_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/030/SRR14580530/SRR14580530_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/030/SRR14580530/SRR14580530_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/031/SRR14580531/SRR14580531_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/031/SRR14580531/SRR14580531_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/032/SRR14580532/SRR14580532_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/032/SRR14580532/SRR14580532_2.fastq.gz

