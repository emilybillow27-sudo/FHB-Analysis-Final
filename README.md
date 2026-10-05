# Genomic Prediction of Fusarium Head Blight Resistance in Wheat

Analysis code for Emily Billow’s master’s thesis in the Colorado State University Wheat Breeding Program.

This research evaluates genomic prediction of Fusarium head blight resistance in hard winter wheat using historical phenotypic and genomic data and CSU field evaluations from 2025 and 2026. It compares historical training data with an expanded training population incorporating eligible 2025 observations to predict performance in a withheld 2026 validation population.

## Traits

- Disease incidence (INC)
- Disease severity (SEV)
- Fusarium-damaged kernels (FDK)
- Deoxynivalenol concentration (DON)

## Repository structure

- `scripts/` — R scripts for phenotypic analysis, genomic prediction, and visualization
- `data/` — Analysis inputs
- `results/` — Generated tables and figures

## Analysis

The analysis includes phenotype adjustment, genomic relationship matrix construction, GBLUP prediction, cross-validation, heritability estimation, and population structure visualization. The repository is under active development; `scripts/BLUEs.R` is the current main pipeline in this checkout.

Run scripts from the repository root. R package requirements and input files are specified in the individual scripts.
