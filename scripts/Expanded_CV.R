# Expanded training population cross-validation
# Emily Billow
# September 2026

# This script loads expanded historical + 2025 multi-environment
# BLUEs prepared by working_script_fixed.R and evaluates genomic
# prediction ability using repeated five-fold cross-validation.

# Run from the project root (FHB Analysis Final).
if (!file.exists("data/FHB_Project_Training_Data.csv")) {
  stop("Set the working directory to the FHB Analysis Final project root.")
}
for (output_dir in c("results/intermediate", "results/tables", "results/figures")) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

library(dplyr)
library(tidyr)
library(sommer)

traits <- c("INC", "SEV", "DON", "FDK")

n_reps <- 100
n_folds <- 5

set.seed(123)

# Load expanded ME BLUEs prepared by working_script_fixed.R
# This script evaluates CV only; it does not reconstruct training phenotypes

expanded_file <- "results/intermediate/blues_me_expanded.rds"
if (!file.exists(expanded_file)) {
  stop("Missing blues_me_expanded.rds. Run the expanded-training preparation section of working_script_fixed.R first.")
}

blues_me_expanded <- readRDS(expanded_file)
GRM <- readRDS("results/intermediate/GRM.rds")

stopifnot(all(c("FullSampleName", "TRAIT", "predicted.value") %in% names(blues_me_expanded)))

# Cross-validation of expanded training population

cv_prediction_list <- list()

for (trait in traits) {
  
  dat_trait <- blues_me_expanded %>%
    dplyr::filter(
      TRAIT == trait,
      !is.na(FullSampleName),
      !is.na(predicted.value)
    ) %>%
    dplyr::transmute(
      FullSampleName = FullSampleName,
      BLUE = predicted.value
    ) %>%
    dplyr::distinct()
  
  dat_trait <- dat_trait %>%
    dplyr::filter(
      FullSampleName %in% rownames(GRM)
    )
  
  GRM_trait <- GRM[
    dat_trait$FullSampleName,
    dat_trait$FullSampleName,
    drop = FALSE
  ]
  
  for (rep_i in seq_len(n_reps)) {
    
    fold_assignment <- sample(
      rep(
        seq_len(n_folds),
        length.out = nrow(dat_trait)
      )
    )
    
    rep_predictions <- list()
    
    for (fold_i in seq_len(n_folds)) {
      
      dat_cv <- dat_trait %>%
        dplyr::mutate(
          FOLD = fold_assignment,
          y_mask = ifelse(
            FOLD == fold_i,
            NA,
            BLUE
          )
        )
      
      dat_cv$FullSampleName <- factor(
        dat_cv$FullSampleName,
        levels = rownames(GRM_trait)
      )
      
      fit_cv <- sommer::mmes(
        fixed = y_mask ~ 1,
        random = ~ vsm(
          ism(FullSampleName),
          Gu = GRM_trait
        ),
        rcov = ~ units,
        data = dat_cv
      )
      
      pred_cv <- predict.mmes(
        fit_cv,
        D = "FullSampleName"
      )$pvals %>%
        as.data.frame()
      
      fold_predictions <- pred_cv %>%
        dplyr::transmute(
          FullSampleName = as.character(
            FullSampleName
          ),
          GBLUP = predicted.value
        ) %>%
        dplyr::inner_join(
          dat_cv %>%
            dplyr::filter(FOLD == fold_i) %>%
            dplyr::transmute(
              FullSampleName = as.character(
                FullSampleName
              ),
              BLUE = BLUE
            ),
          by = "FullSampleName"
        ) %>%
        dplyr::mutate(
          TRAIT = trait,
          REPETITION = rep_i,
          FOLD = fold_i
        )
      
      rep_predictions[[fold_i]] <- fold_predictions
    }
    
    rep_predictions <- dplyr::bind_rows(
      rep_predictions
    )
    
    cv_prediction_list[[length(cv_prediction_list) + 1]] <- rep_predictions
  }
}

cv_predictions_expanded <- dplyr::bind_rows(
  cv_prediction_list
)

# Calculate prediction ability

cv_results_expanded <- cv_predictions_expanded %>%
  dplyr::group_by(
    TRAIT,
    REPETITION
  ) %>%
  dplyr::summarise(
    N = sum(
      complete.cases(BLUE, GBLUP)
    ),
    PREDICTION_ABILITY = cor(
      BLUE,
      GBLUP,
      use = "complete.obs"
    ),
    .groups = "drop"
  )

# Summarize prediction ability across repetitions

cv_summary_expanded <- cv_results_expanded %>%
  dplyr::group_by(TRAIT) %>%
  dplyr::summarise(
    N_REPETITIONS = dplyr::n(),
    MEAN_PREDICTION_ABILITY = mean(
      PREDICTION_ABILITY,
      na.rm = TRUE
    ),
    SD_PREDICTION_ABILITY = sd(
      PREDICTION_ABILITY,
      na.rm = TRUE
    ),
    MIN_PREDICTION_ABILITY = min(
      PREDICTION_ABILITY,
      na.rm = TRUE
    ),
    MAX_PREDICTION_ABILITY = max(
      PREDICTION_ABILITY,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(cv_summary_expanded)

# Cross-validation diagnostics

cv_diagnostics_expanded <- cv_predictions_expanded %>%
  dplyr::group_by(
    TRAIT,
    REPETITION
  ) %>%
  dplyr::summarise(
    N = dplyr::n(),
    N_FOLDS = dplyr::n_distinct(FOLD),
    .groups = "drop"
  )

print(cv_diagnostics_expanded)

cv_fold_sizes_expanded <- cv_predictions_expanded %>%
  dplyr::count(
    TRAIT,
    REPETITION,
    FOLD
  )

print(cv_fold_sizes_expanded)

# Save cross-validation results

saveRDS(
  cv_predictions_expanded,
  "results/intermediate/cv_predictions_expanded.rds"
)


# Pipeline-stage table folders
dir.create("results/tables/cross_validation/expanded", recursive=TRUE, showWarnings=FALSE)
write.csv(
  cv_results_expanded,
  "results/tables/cross_validation/expanded/cv_results_expanded.csv",
  row.names = FALSE
)

write.csv(
  cv_summary_expanded,
  "results/tables/cross_validation/expanded/cv_summary_expanded.csv",
  row.names = FALSE
)