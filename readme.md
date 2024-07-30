## dnam-lung-cancer-pipeline

For nextflow usage Description:

Pipeline for analysing data generated from the DNAm lung cancer screening panel (https://github.com/MRCIEU/dnam-lung-cancer-screening-panel).

### Docker Container Download

To download the required Docker container for the pipeline:

```bash

singularity pull docker://onuroztornaci/dnam_cancer_pipeline:latest

```

### Conda nextflow installation

```bash

conda create -n nextflow_env
conda activate nextflow_env
conda install -c bioconda nextflow

```

### Prepare reference genome 

```bash
GENOME_INDEX=genome/
bash prepare_reference_genome.sh $GENOME_INDEX 

```

[bed_file_path](https://github.com/MRCIEU/dnam-lung-cancer-pipeline/raw/main/data/)

Running time is about 2 hours.

## Usage: 

```
DATA_DIR=data/
PANEL=panel.bed

nextflow main.nf --pipeline bismark \
--data $DATA_DIR \
--genome_folder $GENOME_INDEX \
-resume \
-N [The email address for the report is]

nextflow main.nf --pipeline picard \
--data $DATA_DIR \
--panel $PANEL \
--genome_folder $GENOME_INDEX/hg19.fa \
-resume \
-N [The email address for the report is]
```

*Camda Pipeline*

```
SAMTOOLS_DIR=/bin/

nextflow main.nf \
--data $DATA_DIR \ 
--genome_folder $GENOME_INDEX \
--samtools_path $SAMTOOLS_DIR \
-resume \
-N [The email address for the report is]
```

## To do

* Genome indexing should be a process that runs only if the genome index needs to be created.
```

* Note: It is now possible to have a single directory with all fastq files.
  If new files are generated, just copy them to the directory and
  rerun the pipeline with the '-resume' option
