# DNAm Target Panel Analysis Pipeline

This pipeline enables the analysis of DNA methylation data generated from **any DNAm target panel**.  
While the example provided here uses the [DNAm lung cancer screening panel](https://github.com/MRCIEU/dnam-lung-cancer-screening-panel), the pipeline can be easily adapted for other custom or published DNAm panels by providing an appropriate `panel.csv` file.

---

## 1. Phenotype Matrix (`phenotype.csv`)

- **Sample names** must *exactly* match the corresponding FASTQ file names.  
  *(See the included `phenotype.csv` for an example.)*
- All phenotype variable names should be unique and descriptive.

---

## 2. EWAF Model Matrix (`models.csv`)

### Essentials for Model Matrix Construction

- **Parameters (`p`):** The total number of coefficients estimated by the model, including the intercept and any dummy-coded factors.
- **Minimum complete samples:** You must have at least *(p + 1)* individuals with *no missing values* in all covariates.  
  This is necessary so that the model retains at least one residual degree of freedom *(df = n – rank(X) ≥ 1)*, ensuring valid standard errors and *p*-values.
- **Contrast variables:** For binary indicators derived from multi-category factors, the rows representing the reference group will be set to `NA`.  
  The **effective sample size** for each contrast is the number of rows with non-missing values for that contrast variable.

> **Warning:**  
> If the number of complete samples is less than the number of model parameters, model fitting functions (e.g., `ewaff.sites()` or any ordinary least squares GLM) will return `NA` or `NaN` for t-statistics and p-values due to zero residual degrees of freedom.

**Best practice:**  
Ensure your phenotype file contains *more complete observations* (i.e., with no missing covariate values) than the total number of model parameters.

- The `models.csv` file contains example/default model specifications.
- **Important:**  
  Variable names in your model matrix *must* **exactly match** the variable names in your phenotype file.

---

## 3. DNAm Panel Information (`panel.csv`)

This file specifies the CpG sites included in the DNAm target panel used for your study.  
The format must match the following structure:

| source | details      | chr | start     | end       |
|--------|-------------|-----|-----------|-----------|
| age    | cg20822990  | 1   | 17338766  | 17338766  |
| age    | cg22512670  | 1   | 26855765  | 26855765  |
| age    | cg25410668  | 1   | 28241577  | 28241577  |
| age    | cg04400972  | 1   | 117665053 | 117665053 |
| age    | cg16054275  | 1   | 169556022 | 169556022 |
| age    | cg10501210  | 1   | 207997020 | 207997020 |

**Column descriptions:**
- `source`: Biological context or group (e.g., "age", "smoking") for the CpG site.
- `details`: CpG identifier.
- `chr`: Chromosome number.
- `start`, `end`: Genomic coordinates of the CpG site (usually identical for single-site probes).

> **Note:**  
> To use this pipeline for other panels, simply provide a `panel.csv` with the appropriate CpG list and annotation in the format above.

---

**Example panel:**  
[DNAm lung cancer screening panel](https://github.com/MRCIEU/dnam-lung-cancer-screening-panel) is provided as a template.

**For further details and updates, see the original repository:**  
👉 [https://github.com/MRCIEU/dnam-lung-cancer-screening-panel](https://github.com/MRCIEU/dnam-lung-cancer-screening-panel)
