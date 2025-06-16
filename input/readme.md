# Input File Guide

This guide explains the three files you need: `phenotype.csv`, `models.csv`, and `panel.csv`.

### 1. Your Data File (`phenotype.csv`)

This file contains all the data for each individual in your study.

* The first column **must be `sample_id`**.
* The names in the `sample_id` column **must exactly match** your FASTQ file names.


**Example `phenotype.csv`:**
*(Based on your provided data)*

| sample_id | Sex    | Country.x | Age_Rerc | Cancer_status | BMI_C | Smoke_status | Alc_Lifetime | batch |
| :-------- | :----- | :-------- | :------- | :------------ | :---- | :----------- | :----------- | :---- |
| sample_01 | Female | seven     | 53.9     | Disease       | 25.8  | Former       | 14           | 2     |
| sample_02 | Female | four      | 39.3     | Healthy       | 27.7  | Never        | 2            | 0     |
| sample_03 | Male   | seven     | 47.7     | Disease       | 29.3  | Never        | 15           | 0     |
| sample_04 | Female | four      | 46.1     | Healthy       | 30.4  | Current      | 1            | 0     |

---

### 2. The Model Plan (`models.csv`)

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

### 3. The DNA Panel File (`panel.csv`)

This file lists the CpG sites (DNA locations) your panel measures. This part is standard.

* Your file must have these exact columns: `chr`, `start`, `end`, `source`, `details`.

**Example `panel.csv`:**
| chr  | start    | end      | source | details    |
| :--- | :------- | :------- | :----- | :--------- |
| chr1 | 17338766 | 17338766 | age    | cg20822990 |
| chr1 | 26855765 | 26855765 | age    | cg22512670 |

*To use a different panel, just replace this file with yours in the same format.*

**For an example panel, see:**
👉 [https://github.com/MRCIEU/dnam-lung-cancer-screening-panel](https://github.com/MRCIEU/dnam-lung-cancer-screening-panel)