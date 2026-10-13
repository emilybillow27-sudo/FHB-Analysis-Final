# FHB Analysis

R analysis for Emily Billow's Colorado State University wheat breeding thesis. The project evaluates Fusarium head blight (FHB) traits and compares genomic prediction using historical training data with an expanded population that includes 2025 UIUC observations.

## Traits

- INC: disease incidence (%)
- SEV: disease severity (%)
- FDK: Fusarium-damaged kernels (%)
- DON: deoxynivalenol concentration (ppm)

## Project folders

- `data/`: phenotypic datasets and marker inputs
- `scripts/`: analysis and plotting scripts
- `results/intermediate/`: saved R objects used by subsequent analyses
- `results/tables/`: tabular outputs, organized by analysis
- `results/figures/`: figures, organized by analysis
- `archive/`: older files retained locally and excluded from Git

## Running the analysis

Open `FHB Analysis Final.Rproj` in RStudio and run scripts from the project root so relative paths resolve correctly. Install the R packages required by the scripts and ensure the input files are present in `data/`.

Start with the main script to prepare adjusted genotype estimates, construct the genomic relationship matrix, and generate forward predictions:

```r
source("scripts/working_script_fixed.R")
```

Then run the analyses needed for the current task:

```r
source("scripts/Historical_CV.R")
source("scripts/Expanded_CV.R")
source("scripts/compare_training_scenarios_2026.R")
source("scripts/Heritabilities.R")
source("scripts/PCA.R")
source("scripts/Distributions.R")
```

Other plotting scripts summarize cross-validation, trait correlations, population comparisons, and observed versus predicted values. If marker inputs, phenotype inputs, or the main analysis change, rerun dependent analyses before using their outputs.

## Current data status

The 2026 results are preliminary. UIUC and SDSU DON measurements are pending, and SDSU FDK is currently unavailable. Missing measurements remain missing and must not be interpreted as zero. Current 2026 DON predictions are evaluated against UNL observations only.

## Interpretation

Prediction ability is the Pearson correlation between predicted genomic breeding values and observed genotype estimates. Historical and expanded cross-validation evaluate different populations. Forward prediction compares the training scenarios using the same validation lines, with their phenotypic records excluded from training.
