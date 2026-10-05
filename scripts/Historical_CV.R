# Historical population genomic prediction cross-validation
# Author: Emily Billow
# Purpose: Evaluate GBLUP prediction ability within the historical FHB
# training population using 100 repetitions of 5-fold cross-validation.

# Run from the project root (FHB Analysis Final).
if (!file.exists("data/FHB_Project_Training_Data.csv")) {
  stop("Set the working directory to the FHB Analysis Final project root.")
}
for (output_dir in c("results/intermediate", "results/tables", "results/figures")) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

library(dplyr)
library(sommer)

traits <- c("INC", "SEV", "DON", "FDK")

n_reps <- 100
n_folds <- 5

set.seed(123)

# Load finalized historical multi-environment BLUEs and genomic relationship matrix
blues_me_train <- readRDS("results/intermediate/blues_me_train.rds")
GRM <- readRDS("results/intermediate/GRM.rds")

# Store out-of-fold predictions from each trait and repetition
cv_prediction_list <- list()

for (trait in traits) {
  
  message("Historical CV for trait: ", trait)
  
  # Select historical phenotypes for the current trait
  dat_trait <- blues_me_train %>%
    dplyr::filter(
      TRAIT == trait,
      !is.na(FullSampleName),
      !is.na(predicted.value)
    ) %>%
    dplyr::select(
      FullSampleName,
      BLUE = predicted.value
    ) %>%
    dplyr::distinct(FullSampleName, .keep_all = TRUE)
  
  # Retain only genotypes represented in the genomic relationship matrix
  dat_trait <- dat_trait %>%
    dplyr::filter(FullSampleName %in% rownames(GRM))
  
  message("Number of historical genotypes: ", nrow(dat_trait))
  
  # Subset the genomic relationship matrix to the current population
  GRM_trait <- GRM[
    dat_trait$FullSampleName,
    dat_trait$FullSampleName,
    drop = FALSE
  ]
  
  for (rep_i in seq_len(n_reps)) {
    
    message(
      "Trait: ", trait,
      " | repetition: ", rep_i, "/", n_reps
    )
    
    # Randomly assign genotypes to five approximately equal folds
    fold_assignment <- sample(
      rep(
        seq_len(n_folds),
        length.out = nrow(dat_trait)
      )
    )
    
    rep_predictions <- list()
    
    for (fold_i in seq_len(n_folds)) {
      
      # Mask phenotypes for genotypes in the validation fold
      dat_cv <- dat_trait %>%
        dplyr::mutate(
          FOLD = fold_assignment,
          y_mask = ifelse(
            FOLD == fold_i,
            NA_real_,
            BLUE
          )
        )
      
      # Match genotype factor levels to the genomic relationship matrix
      dat_cv$FullSampleName <- factor(
        dat_cv$FullSampleName,
        levels = rownames(GRM_trait)
      )
      
      # Fit GBLUP using phenotypes from the other four folds
      fit_cv <- sommer::mmes(
        fixed = y_mask ~ 1,
        random = ~ sommer::vsm(
          sommer::ism(FullSampleName),
          Gu = GRM_trait
        ),
        rcov = ~ units,
        data = dat_cv,
        dateWarning = FALSE
      )
      
      # Obtain genomic predictions
      pred_cv <- sommer::predict.mmes(
        fit_cv,
        D = "FullSampleName"
      )$pvals %>%
        dplyr::select(
          FullSampleName,
          GBLUP = predicted.value
        )
      
      # Retain predictions for the withheld validation fold
      fold_predictions <- dat_cv %>%
        dplyr::filter(FOLD == fold_i) %>%
        dplyr::select(
          FullSampleName,
          BLUE,
          FOLD
        ) %>%
        dplyr::mutate(
          FullSampleName = as.character(FullSampleName)
        ) %>%
        dplyr::left_join(
          pred_cv %>%
            dplyr::mutate(
              FullSampleName = as.character(FullSampleName)
            ),
          by = "FullSampleName"
        ) %>%
        dplyr::mutate(
          TRAIT = trait,
          REPETITION = rep_i
        )
      
      rep_predictions[[fold_i]] <- fold_predictions
    }
    
    # Combine out-of-fold predictions from all five folds
    rep_predictions <- dplyr::bind_rows(rep_predictions)
    
    cv_prediction_list[[length(cv_prediction_list) + 1]] <- rep_predictions
  }
}

# Combine predictions across all traits and repetitions
cv_predictions_hist <- dplyr::bind_rows(cv_prediction_list)

# Calculate prediction ability for each repetition
cv_results_hist <- cv_predictions_hist %>%
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

# Summarize prediction ability across the 100 repetitions
cv_summary_hist <- cv_results_hist %>%
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

# Check that 100 prediction-ability estimates were generated for each trait
cv_results_hist %>%
  dplyr::count(TRAIT)

# Check the number of genotypes predicted in each repetition
cv_results_hist %>%
  dplyr::select(
    TRAIT,
    REPETITION,
    N
  )

# Check fold sizes
cv_predictions_hist %>%
  dplyr::count(
    TRAIT,
    REPETITION,
    FOLD
  )

# View final summary
print(cv_summary_hist)

# Save genotype-level out-of-fold predictions
saveRDS(
  cv_predictions_hist,
  "results/intermediate/historical_cv_predictions.rds"
)

# Save prediction ability for each repetition
write.csv(
  cv_results_hist,
  "results/tables/historical_cv_results.csv",
  row.names = FALSE
)

# Save trait-level summary
write.csv(
  cv_summary_hist,
  "results/tables/historical_cv_summary.csv",
  row.names = FALSE
)