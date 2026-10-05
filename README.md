# Genomic Prediction of Fusarium Head Blight Resistance in Wheat

Emily Billow’s master’s thesis analysis in the Colorado State University Wheat Breeding Program. The study evaluates historical and expanded training populations for prediction of incidence, severity, Fusarium-damaged kernels, and DON in CSU wheat material evaluated in 2025 and 2026.

## Scripts

- `scripts/working_script_fixed.R` — Main pipeline: phenotype estimates, GRM, and forward prediction
- `scripts/Historical_CV.R` — Historical training population cross-validation
- `scripts/Expanded_CV.R` — Expanded training estimates and cross-validation
- `scripts/PCA.R` — Population structure and scree plots
- `scripts/scatterplots.R` — Observed versus predicted plots

Run scripts from the project root. The main pipeline reads the saved expanded-training estimates; refresh those with `Expanded_CV.R` after updating historical or 2025 inputs, then rerun the main pipeline for the corresponding 2026 comparison. CV scripts retain 100 repetitions of five-fold validation. PCA and scatterplots read saved inputs.

## Workspace

- `data/` — Original analysis inputs
- `results/intermediate/` — Saved phenotype estimates, marker inputs, GRM, and CV prediction objects
- `results/tables/` — Prediction summaries, PCA exports, and retained research summary tables
- `results/figures/` — PCA and prediction figures from the active scripts

Superseded scripts and alternate outputs are preserved locally under `archive/2026-10-05-cleanup/` and remain recoverable in Git history. Existing result files were organized without recalculating the statistical analyses. Large binary inputs and saved R objects use Git LFS.
