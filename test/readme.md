## Minimal Test Dataset

- Includes a **small, reproducible dataset** derived from real WGBS data ([ENA PRJNA730913](https://www.ebi.ac.uk/ena/browser/view/PRJNA730913)), designed for fast prototyping and demonstration.
- All files are based on a subset of the original study, selected to minimize runtime and disk usage.
- Every file can be **fully regenerated** using the provided scripts.  
- No private or protected data is included — everything is public and open-access.

---

## Key contents

- [input/](input): files for test dataset analysis
    - [input/panel.csv](input/panel.csv): panel of capture regions for the test dataset (created by scripts/select-test-panel.sh)
    - [input/models.csv](input/models.csv): models for testing associations in the test dataset
    - [input/phenotypes.csv](input/phenotypes.csv): phenotype data for testing associations in the test dataset
- [data/raw](data/raw): fastq files for the test dataset (created by scripts/create-test-fastq-files.sh)
- [data/blood-cell-type-reference.csv](data/blood-cell-type-reference.csv): blood cell type DNA methylation reference dataset
- [data/genome-reference/](data/genome-reference): genome reference for aligning the test dataset (created by scripts/create-test-reference.sh) 
- [scripts/create-test-dataset.sh](scripts/create-test-dataset.sh): script to recreate the test dataset

---

## Recreate the test dataset

```
bash scripts/create-test-dataset.sh 
```

Scripts that depend on tools that are not part of a
typical linux environment use apptainer containers.
Consequently, the only dependency of this pipeline is apptainer.


