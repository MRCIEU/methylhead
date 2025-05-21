# Information about Phenotype and EWAF Model Matrices

## Phenotype Matrix (`phenotype.csv`)

- **Sample names** must exactly match the corresponding FASTQ file names.
    - *Example: see `phenotype.csv`.*
- All phenotype variable names should be unique and descriptive.

## EWAF Model Matrix (`models.csv`)

### Key Points for Building Your Model Matrix

- **Parameters (p):** Total number of coefficients estimated by the model (including the intercept and any dummy-coded factors).
- **Minimum complete samples:** The minimum number of individuals required with *no missing values* in all covariates, so that the model has at least one residual degree of freedom *(i.e., df = n – rank(X) ≥ 1)*. This ensures that valid standard errors and *p*-values can be computed*.
- **Contrast variables:** When using binary indicators derived from multi-category factors, the rows for the reference group will be set to `NA`. 
- **The *effective sample size*** for these contrasts is the number of rows with non-missing values for the respective contrast variable.

> **Warning:**  
> If you have fewer complete samples than the number of model parameters, functions like `ewaff.sites()` (or any ordinary least squares GLM) will return `NA` or `NaN` for t-statistics and p-values because there are no residual degrees of freedom.

**Therefore:**  
Make sure your phenotype file contains more *complete* observations (with no missing covariate values) than the total number of model parameters.

- The `models.csv` file contains example/default model specifications.
- **Important:** Variable names in the model matrix *must* exactly match the variable names in your phenotype file.
