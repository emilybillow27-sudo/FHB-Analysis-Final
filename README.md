# Genomic Prediction of Fusarium Head Blight Resistance in Wheat

Emily Billow’s master’s thesis analysis in the Colorado State University Wheat Breeding Program. The study compares historical and expanded training populations for prediction of incidence, severity, Fusarium-damaged kernels, and DON in CSU wheat material evaluated in 2025 and 2026.

## Run the analysis

Install Git LFS before cloning and run `git lfs pull` to retrieve data. Restore the recorded R dependencies with `renv::restore()`, then run:

```sh
Rscript run_analysis.R --stage all
```

`config.R` defines inputs, traits, marker thresholds, random seed, and the default 100 repetitions of five-fold cross-validation. Individual stages can be run with `--stage prepare`, `prediction`, `cv`, `summaries`, or `figures`. A shorter validation run uses `--repetitions 2`; it is not a final analysis.

The pipeline prepares harmonized phenotype estimates and marker data, constructs both training populations, evaluates forward and cross-validation predictions, and generates summaries and figures. All stages read explicit saved inputs rather than relying on an existing R session.

## Files and outputs

- `data/` — Phenotype and marker inputs
- `scripts/` — Shared functions and preparation, prediction, summary, and figure modules
- `results/pipeline/` — New analysis outputs, model diagnostics, and run manifests; generated locally
- `results/` — Previously committed results retained for reference

Use `--output path/to/results` for a separate run. Missing location-specific traits are skipped until measurements are available. Heritability tables flag singular fits; these require interpretation before reporting.
