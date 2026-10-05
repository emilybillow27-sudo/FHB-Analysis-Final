# FHB Analysis Final

Emily Billow’s master’s thesis analysis of genomic prediction for Fusarium head blight resistance in hard winter wheat at Colorado State University.

The project compares historical training data with historical data expanded using eligible 2025 CSU phenotypes, then evaluates prediction of a withheld 2026 validation population. Traits are incidence (INC), severity (SEV), Fusarium-damaged kernels (FDK), and deoxynivalenol (DON).

## Repository synchronization

The previous GitHub README identified `scripts/BLUEs.R` as the main in-progress pipeline. The current Desktop working copy contains newer scripts, including `scripts/working_script_fixed.R`, both root CV scripts, and `working_heritabilities.R`, that have not yet been committed to GitHub. The inventory below describes that local working copy. This documentation update does not upload those scripts, data, or results. Synchronizing and validating them is an open task.

## Project navigation

- [Project tasks](PROJECT_TASKS.md): priorities, milestones, and completion checks.
- [Analysis workflow](ANALYSIS_WORKFLOW.md): scripts, dependencies, and reproducibility gaps.
- `data/`: local phenotype and marker inputs. Availability and permission to share must be established before uploading new research data.
- `scripts/`: analysis and plotting scripts.
- `results/`: existing generated estimates, prediction summaries, tables, and figures. These are provisional until input provenance and model checks are confirmed.

## Current local working scripts

| Script | Purpose |
| --- | --- |
| `scripts/working_script_fixed.R` | Phenotypic estimates, GRM, and 2025/2026 forward prediction |
| `Historical_CV.R` | Historical five-fold CV repeated 100 times |
| `Expanded_CV.R` | Expanded phenotype estimates and five-fold CV repeated 100 times |
| `scripts/PCA.R` | PCA and scree plots |
| `working_heritabilities.R` | Historical, 2025, and UNL/UIUC 2026 heritability estimates |
| `scripts/Summary Tables.R` | Historical, 2025, and UNL/UIUC 2026 phenotype summaries |
| `scatterplots.R` | Prediction plotting script; confirm current inputs before use |
| `plot_trait_correlations_slide.R` | Presentation trait-correlation plots |

Several overlapping older scripts remain in `scripts/`. Select the authoritative version for each analysis before archiving alternatives; filenames alone do not establish which version is final.

## Running the analysis

Use the repository root as the R working directory. The current scripts are not yet a clean end-to-end pipeline: the main script needs an expanded-training output that is itself built from main-script outputs. See [Analysis workflow](ANALYSIS_WORKFLOW.md) before execution. Do not source all scripts blindly.

R packages directly referenced by the reviewed core scripts include `dplyr`, `tidyr`, `purrr`, `sommer`, `gaston`, `stringr`, `tibble`, `ggplot2`, `lme4`, `readr`, and `e1071`. Package versions have not yet been recorded or locked.

## GitHub organization

Use the repository task checklist for the working plan. Convert individual unchecked items into GitHub issues when remote access is available. Suggested labels are `analysis`, `data`, `writing`, `presentation`, `reproducibility`, and `blocked`.

Commit code and documentation in small, descriptive changes. Review raw data, large marker files, generated binary results, and existing local deletions separately before staging. Do not commit R session files or credentials.
