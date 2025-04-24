# Epigenetic and Lifestyle Correlates of a Novel Healthy Ageing Metric

This project analyzes LBC1936 Wave 1 data to explore how epigenetic clocks and lifestyle factors relate to a latent "Healthy Ageing" (HA) metric, constructed from cognitive, physical, biological, and quality of life (QOL) indicators.

## Workflow Overview

### `01_Data_Prep/`
- **QC of selected variables**: histograms, Shapiro-Wilk tests
- **Outlier handling**: recode outliers to NA
- **Summary statistics**

### `02_Regressions/`
- **Regression adjustment**: remove effects of age, sex (and height) from variables

### `03_Plots/`
- **Further QC**: bivariate and density plots to detect missed outliers

### `04_SEM/`
- **Latent variable modeling**:
  - Domains: cognitive, physical, biological, QOL
  - Higher-order latent variable: Healthy Ageing (HA)
- **Outputs**:
  - Path diagrams
  - Scree plots
  - Correlation matrices
  - Extracted latent scores

### `05_Linear_Models/`
- **Modeling**:
  - Linear models predicting domains and HA
  - Epigenetic clocks and lifestyle factors as predictors
- **QC**: of clocks data
- **Diagnostics**: posterior predictive check, linearity, homogeneity of variance, normality of residuals
- **Reporting**: forest plots for model results
  

## Notes
- Uses `LBC1936` data — not included in this repo due to data sensitivity
- Paths are relative and assume project is run from the wd `/Volumes/igmm/marioni-lab/Ella_Shuttleworth`
- This project requires several R packages, which will be installed automatically when you run the scripts.
