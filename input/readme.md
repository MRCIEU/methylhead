# Input file guide

This guide explains the four main input files required by the pipeline: `phenotype.csv`, `models.csv`, `panel.csv` and `samplesheet.csv`.

### 1. Sample phenotype information (`phenotype.csv`)

This file contains all the data for each individual in your study.

* The first column **must be `sample_id`**.


**Example `phenotype.csv`:**

| sample_id | Sex    | Country.x | Age_Rerc | Cancer_status | BMI_C | Smoke_status | Alc_Lifetime | batch |
| :-------- | :----- | :-------- | :------- | :------------ | :---- | :----------- | :----------- | :---- |
| sample_01 | Female | seven     | 53.9     | Disease       | 25.8  | Former       | 14           | 2     |
| sample_02 | Female | four      | 39.3     | Healthy       | 27.7  | Never        | 2            | 0     |
| sample_03 | Male   | seven     | 47.7     | Disease       | 29.3  | Never        | 15           | 0     |
| sample_04 | Female | four      | 46.1     | Healthy       | 30.4  | Current      | 1            | 0     |

---

### 2. Models for association testing (`models.csv`)

This file is your "recipe book" for all the statistical analyses you want to run. Each row is one model.

* **Most Important:** All variable names used in the `model` column (like `Alc_Lifetime`, `Age_Rerc`, `Sex`, etc.) **must exist** as column headers in your `phenotype.csv` file.

**Explanation of Columns:**
* `name`: A short, unique name for each analysis.
* `var`: The main variable you are interested in for that analysis.
* `model`: The exact statistical formula to run. It usually follows the format `methylation ~ your_main_variable + variable_to_adjust_for_1 + variable_to_adjust_for_2`.

**Example `models.csv`:**
*(Based on your provided models)*

| name       | var            | model                                                     |
| :--------- | :------------- | :-------------------------------------------------------- |
| alcohol    | Alc_Lifetime   | methylation~Alc_Lifetime+Smoke_status+BMI_C+Age_Rerc+Sex  |
| bmi        | BMI_C          | methylation~BMI_C+Smoke_status+Age_Rerc+Sex               |
| cancer     | Cancer_status  | methylation~Cancer_status+Age_Rerc                        |
| cancer_adj | Cancer_status  | methylation~Cancer_status+Age_Rerc+Sex+Smoke_status+Alc_Lifetime |

---

### 3. Genomic targets (`panel.csv`)

This file provides the genomic coordinates of targeted regions. 
These are used by the pipeline to evaluate the targeting assay.

* Your file must have these exact columns: `chr`, `start`, `end`, `source`, `details`.

**Example `panel.csv`:**
| chr  | start    | end      | source | details    |
| :--- | :------- | :------- | :----- | :--------- |
| chr1 | 17338766 | 17338766 | age    | cg20822990 |
| chr1 | 26855765 | 26855765 | age    | cg22512670 |

### 4. Samplesheet (`samplesheet.csv`)

This file lists the paired fastq input files for each sample.

**Explanation of Columns:**
* `sample_id`: Sample identifier, must match those in `phenotype.csv`.
* `read1`: Left end reads
* `read2`: Right end reads 

**Example `samplesheet.csv`:**

|sample_id|read1|read2|
| :--- | :------- | :------- | 
|sample_01|/data/methyl-seq/raw/SRR14580504_L001_1.fastq.gz|/data/methyl-seq/raw/SRR14580504_L001_2.fastq.gz|
|sample_01|/data/methyl-seq/raw/SRR14580504_L002_1.fastq.gz|/data/methyl-seq/raw/SRR14580504_L002_2.fastq.gz|
|sample_02|/data/methyl-seq/raw/SRR14580505_L001_1.fastq.gz|/data/methyl-seq/raw/SRR14580505_L001_2.fastq.gz|
|sample_02|/data/methyl-seq/raw/SRR14580505_L002_1.fastq.gz|/data/methyl-seq/raw/SRR14580505_L002_2.fastq.gz|

