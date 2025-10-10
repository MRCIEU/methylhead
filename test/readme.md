## Minimal Test Dataset

The test dataset was derived from a real methyl-seq data 
([ENA PRJNA730913](https://www.ebi.ac.uk/ena/browser/view/PRJNA730913))
using scripts in the [scripts/](scripts) folder.

Processing this dataset using the methylhead pipeline takes about 10 minutes.

---

## Key contents

- [data/panel.csv](input/panel.csv): panel of capture regions for the test dataset (created by scripts/select-test-panel.sh)
- [data/models.csv](input/models.csv): models for testing associations in the test dataset
- [data/phenotypes.csv](input/phenotypes.csv): phenotype data for testing associations in the test dataset
- [data/raw](data/raw): fastq files for the test dataset (created by scripts/create-test-fastq-files.sh)
- [data/samplesheet.csv](data/samplesheet.csv): list of all fastq files with corresponding sample identifiers
- [data/blood-cell-type-reference.csv](data/blood-cell-type-reference.csv): blood cell type DNA methylation reference dataset
- [scripts/create-test-dataset.sh](scripts/create-test-dataset.sh): script to recreate the test dataset
- [nextflow.config](nextflow.config): nextflow configuration file

---

## Running the pipeline

Please first edit [nextflow.config](nextflow.config) to specify 
the location of your indexed hg19 genome (see line 6, 'genome_dir'). 
See pipeline [instructions](../readme.md) for downloading and indexing the genome.

```
nextflow -C nextflow.config run ../main.nf -resume
```

As defined by `nextflow.config`, 
pipeline outputs will be appear locally in the 'output' folder, 
and the nextflow work outputs in the 'work' folder.

---

## Recreate the test dataset

*This step has already been performed and all outputs are contained the [data/](data) folder.*

```
bash scripts/create-test-dataset.sh PATH/TO/ORIGINAL/DATA PATH/TO/HG19/GENOME/INDEX
```

* `PATH/TO/ORIGINAL/DATA` The location where the original dataset should
be (or has been) downloaded and aligned.
* `PATH/TO/HG19/GENOME/INDEX` The location where the human genome (hg19) should be (or has been) downloaded and indexed.



