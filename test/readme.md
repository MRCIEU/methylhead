## Minimal Test Dataset

- Includes a **small, reproducible dataset** derived from real WGBS data ([ENA PRJNA730913](https://www.ebi.ac.uk/ena/browser/view/PRJNA730913)), designed for fast prototyping and demonstration.
- All files are based on a subset of the original study, selected to minimize runtime and disk usage.
- Every file can be **fully regenerated** using the provided scripts.  
- No private or protected data is included — everything is public and open-access.

---
### Contents

- [input/](input): files for test dataset analysis
    - [input/panel.csv](input/panel.csv): panel of capture regions for the test dataset (created by src/select-test-panel.sh)
    - [input/models.csv](input/models.csv): models for testing associations in the test dataset
    - [input/phenotypes.csv](input/phenotypes.csv): phenotype data for testing associations in the test dataset
- [reference/](reference): genome reference for aligning the test dataset (created by src/create-test-reference.sh)
- [dataset/](dataset): fastq files for the test dataset (created by src/create-test-fastq-files.sh)
- [src/](src): scripts to recreate the test dataset
    - [src/create-test-dataset.sh](src/create-test-dataset.sh): script to recreate the test dataset, uses the scripts listed below
    - [src/download-dataset.sh](src/download-dataset.sh): downloads the dataset
    - [src/create-reference.sh](src/create-reference.sh): creates reference for the full genome
	- [src/align-dataset.sh](src/align-dataset.sh): aligns reads from the original dataset to the genome (hg19) with BWAmeth 
    - [src/select-test-panel.sh](src/select-test-panel.sh): creates panel of regions for test dataset
    - [src/create-test-fastq-files.sh](src/create-test-fastq-files.sh): creates fastq files for test dataset
    - [src/create-test-reference.sh](src/create-test-reference.sh): creates reference for the test genome

---

## Recreate the test dataset

```
bash src/create-test-dataset.sh
```

**Requirements**

- **BWAmeth**
- **BEDtools**
- **Samtools**
- **Python 3**
- Standard UNIX tools: `bash`, `awk`


