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

```
GENOMES=[path to converted genome index]

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
mv hg19.fa.gz ${GENOMES}
bismark_genome_preparation --path_to_aligner /data/ms13525/lib/miniconda2/envs/Bismark/bin  --verbose ${GENOMES}
wget https://www.twistbioscience.com/sites/default/files/resources/2022-06/covered_targets_Twist_Methylome_hg19_annotated_collapsed_final.bed.zip
mv covered_targets_Twist_Methylome_hg19_annotated_collapsed_final.bed.zip ${GENOMES}
samtools faidx hg19.fa
picard CreateSequenceDictionary REFERENCE=hg19.fa OUTPUT=hg19.fa.dict
picard BedToIntervalList \
I=covered_targets_Twist_Methylome_hg19_annotated_collapsed_final.bed \
O=covered_targets_Twist_Methylome_hg19_annotated_collapsed_final \
SD=hg19.fa.dict

```
Running time is about 2 hours.

## Usage: 

```
nextflow main.nf --pipeline bismark \
--data "[fastq path]" \
--genome_folder [genome index path] \
-resume \
-N [The email address for the report is]

nextflow main.nf --pipeline picard \
--data "[fastq path]" 
--intervals covered_targets_Twist_Methylome_hg19_annotated_collapsed_final \
--genome_folder [BWA genome index path] \
-resume \
-N [The email address for the report is]
```

*Camda Pipeline*

```
nextflow main.nf --data "[fastq path]" \ 
--genome_folder [genome folder path] \
--samtools_path [samtools fodler path] \
-resume \
-N [The email address for the report is]
```

## To do

* Genome indexing should be a process that runs only if the genome index needs to be created.
```

* Note: It is now possible to have a single directory with all fastq files.
  If new files are generated, just copy them to the directory and
  rerun the pipeline with the '-resume' option
