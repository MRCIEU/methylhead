## Small Example Dataset

This dataset was derived from a real methyl-seq data 
([ENA PRJNA730913](https://www.ebi.ac.uk/ena/browser/view/PRJNA730913)).

Processing this dataset takes about 10 minutes.

---

## Key contents

- [data/panel.csv](input/panel.csv): panel of capture regions for the dataset (created by scripts/select-panel.sh)
- [data/models.csv](input/models.csv): models for testing associations in the example dataset
- [data/phenotypes.csv](input/phenotypes.csv): phenotype data for testing associations in the example dataset
- [data/raw](data/raw): fastq files for the example dataset (created by scripts/create-fastq-files.sh)
- [data/samplesheet.csv](data/samplesheet.csv): list of all fastq files with corresponding sample identifiers
- [data/blood-cell-type-reference.csv](data/blood-cell-type-reference.csv): blood cell type DNA methylation reference dataset
- [scripts/create-dataset.sh](scripts/create-dataset.sh): script to recreate the example dataset
- [config.yml](config.yml): snakemake configuration file

---

## Running the pipeline


Please first edit [config.yml](config.yml) to specify 
the location of your indexed hg19 genome ('paths:genome'). 

The entire pipeline can be applied to the dataset 
using the following command. 

```
bash analyse-dataset.sh
```

---

## Recreate the dataset

*The dataset has already been created, so there is no need to repeat it*.
*All outputs are contained the [data/](data) folder*.

```
bash scripts/create-dataset.sh PATH/TO/ORIGINAL/DATA PATH/TO/GENOME/INDEX
```

* `PATH/TO/ORIGINAL/DATA` The location where the original dataset should
be (or has been) downloaded and aligned.
* `PATH/TO/GENOME/INDEX` The location where the human genome (hg19) should be 
(or has been) downloaded and indexed.



