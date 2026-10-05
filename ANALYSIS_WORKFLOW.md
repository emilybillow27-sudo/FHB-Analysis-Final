# Analysis workflow

## Status of this review

This inventory is based on reading local scripts and listing existing outputs on October 5, 2026. Models were not rerun, saved outputs were not validated, and GitHub synchronization was not verified. Existing results do not prove that the latest script version produced them.

## Dependency map

1. `scripts/working_script_fixed.R` reads historical and 2025 phenotype CSVs and creates `blues_se_train.rds`, `blues_se_test.rds`, and `blues_me_train.rds`.
2. It reads `data/fhb_analysis_2026_production_final.vcf.gz`, applies MAF >= 0.05 and heterozygosity <= 0.10 filters, creates `GRM.rds`, harmonizes in-memory marker IDs, and calculates 2025 forward predictions.
3. It reads the three 2026 location CSVs, identifies 60 replicated validation genotypes, and creates single- and across-location estimates.
4. `Expanded_CV.R` reads the historical/2025 single-environment RDS files and GRM, excludes validation names, builds `blues_me_expanded.rds`, and performs expanded CV.
5. The final portion of the main script reads `blues_me_expanded.rds` and compares historical versus expanded 2026 forward prediction.
6. `Historical_CV.R` reads `blues_me_train.rds` and `GRM.rds` independently.
7. `scripts/PCA.R` expects in-memory `geno_mat`, `blues_me_train`, and `blues_se_test`. The reviewed main script does not define `geno_mat`.
8. Heritability and summary scripts read phenotype CSVs separately. Their data selection must be reconciled with the core pipeline.

This map documents dependencies; it is not an executable run order. Resolve the cycle before introducing a single runner.

## Recommended pipeline stages

- Prepare inputs and verify phenotype/marker identifiers.
- Calculate phenotypic estimates and GRM; save harmonized inputs.
- Construct expanded training estimates with explicit validation exclusions.
- Perform historical and expanded CV.
- Evaluate 2025 and 2026 forward predictions, including location-specific evaluations.
- Generate heritability, phenotype summaries, PCA, and trait relationships.
- Export final thesis/presentation tables and figures with a run manifest.

## Checks before finalizing results

- Save phenotype files after identifier harmonization, including the historical single-environment estimates consumed by expanded training.
- Check uniqueness of genotype-to-marker mappings and prevent duplicate CV rows per marker ID.
- Enforce validation exclusions using both biological genotype and canonical marker identifiers.
- Confirm all intended validation IDs occur in the GRM; an ID count alone does not verify matching.
- Record actual usable sample counts by trait, population, and location.
- Confirm historical and expanded forward comparisons use identical observed/predicted validation pairs.
- Review convergence, variance components, missing predictions, and non-finite/zero standard errors before inverse-variance weighting.
- Review phenotype-adjustment choices and validation design with Zach, including preprocessing before CV.
- Define PCA marker coding and imputation explicitly.
- Reconcile heritability scripts and summaries: the reviewed root heritability script and summary script omit SDSU; the root heritability script limits 2026 traits to INC, SEV, and FDK.
- Record software versions, input versions, seed, run date, and output provenance.

## Main output groups

- Phenotype estimates: `results/blues_*.rds`
- Relationship matrix: `results/GRM.rds`
- Historical CV: `results/historical_cv_*.csv` and prediction RDS
- Expanded CV: `results/cv_results_expanded.csv`, `cv_summary_expanded.csv`, and prediction RDS
- Forward prediction: `results/forward_predictions_2025.csv`, `forward_predictions_2026.csv`, and associated summaries/comparison
- Supporting summaries and plots: `results/tables/` and `results/figures/`

Select final outputs explicitly rather than treating every similarly named file as current.
