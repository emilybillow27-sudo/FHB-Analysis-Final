# Expanded training population cross-validation
# Emily Billow
# September 2026

# This script calculates multi-environment BLUEs for the expanded
# historical + 2025 training population and evaluates genomic
# prediction ability using repeated five-fold cross-validation.

library(dplyr)
library(tidyr)
library(sommer)

traits <- c("INC", "SEV", "DON", "FDK")

n_reps <- 100
n_folds <- 5

set.seed(123)

# Load data

blues_se_train <- readRDS(
  "results/blues_se_train.rds"
)

blues_se_test <- readRDS(
  "results/blues_se_test.rds"
)

GRM <- readRDS(
  "results/GRM.rds"
)

unl_2026_raw <- read.csv(
  "data/FHB_Project_Testing_Data_UNL_2026.csv",
  stringsAsFactors = FALSE
)

# Identify 2026 validation genotypes

selected_2026 <- unl_2026_raw %>%
  dplyr::filter(CATEGORY == "Duplicate") %>%
  dplyr::distinct(ID) %>%
  dplyr::pull(ID)

stopifnot(length(selected_2026) == 60)

print(length(selected_2026))

# Prepare historical single-environment BLUEs

hist_se_expanded <- blues_se_train %>%
  dplyr::filter(
    !germplasmName %in% selected_2026
  ) %>%
  dplyr::transmute(
    germplasmName = germplasmName,
    FullSampleName = FullSampleName,
    TRAIT = TRAIT,
    environment = studyName,
    predicted.value = predicted.value,
    std.error = std.error,
    adjusted = adjusted,
    source = "Historical"
  )

# Prepare 2025 single-environment BLUEs

test_2025_se_expanded <- blues_se_test %>%
  dplyr::filter(
    !ID %in% selected_2026
  ) %>%
  dplyr::transmute(
    germplasmName = ID,
    FullSampleName = FullSampleName,
    TRAIT = TRAIT,
    environment = SUB_NURNAME,
    predicted.value = predicted.value,
    std.error = std.error,
    adjusted = adjusted,
    source = "2025"
  )

# Combine historical and 2025 single-environment BLUEs

blues_se_expanded <- dplyr::bind_rows(
  hist_se_expanded,
  test_2025_se_expanded
)

# Confirm 2026 validation genotypes are excluded

selected_remaining <- intersect(
  selected_2026,
  unique(blues_se_expanded$germplasmName)
)

stopifnot(length(selected_remaining) == 0)

# Summarize expanded training population

expanded_population_summary <- blues_se_expanded %>%
  dplyr::group_by(source, TRAIT) %>%
  dplyr::summarise(
    N_RECORDS = dplyr::n(),
    N_GENOTYPES = dplyr::n_distinct(germplasmName),
    N_ENVIRONMENTS = dplyr::n_distinct(environment),
    .groups = "drop"
  )

print(expanded_population_summary)

# Calculate precision-weighted multi-environment BLUEs

expanded_blue_list <- list()

for (trait in traits) {
  
  dat_trait <- blues_se_expanded %>%
    dplyr::filter(
      TRAIT == trait,
      !is.na(predicted.value),
      !is.na(germplasmName),
      !is.na(FullSampleName)
    )
  
  genotype_names <- dat_trait %>%
    dplyr::select(
      germplasmName,
      FullSampleName
    ) %>%
    dplyr::distinct()
  
  dat_model <- dat_trait %>%
    dplyr::select(
      germplasmName,
      environment,
      predicted.value,
      std.error
    ) %>%
    dplyr::rename(
      y = predicted.value,
      se = std.error
    )
  
  dat_model$germplasmName <- factor(
    dat_model$germplasmName
  )
  
  dat_model$environment <- factor(
    dat_model$environment
  )
  
  if (all(is.na(dat_model$se))) {
    
    mod <- sommer::mmes(
      fixed = y ~ germplasmName,
      random = ~ environment,
      rcov = ~ units,
      data = dat_model
    )
    
  } else {
    
    mean_se <- mean(
      dat_model$se,
      na.rm = TRUE
    )
    
    dat_model <- dat_model %>%
      dplyr::mutate(
        se_weight = ifelse(
          is.na(se),
          mean_se,
          se
        )
      )
    
    W <- diag(
      1 / dat_model$se_weight^2
    )
    
    mod <- sommer::mmes(
      fixed = y ~ germplasmName,
      random = ~ environment,
      rcov = ~ units,
      W = W,
      data = dat_model
    )
  }
  
  pred <- predict.mmes(
    mod,
    D = "germplasmName"
  )$pvals %>%
    as.data.frame()
  
  pred <- pred %>%
    dplyr::transmute(
      germplasmName = as.character(germplasmName),
      predicted.value = predicted.value,
      std.error = std.error
    ) %>%
    dplyr::left_join(
      genotype_names,
      by = "germplasmName"
    ) %>%
    dplyr::mutate(
      TRAIT = trait
    )
  
  expanded_blue_list[[trait]] <- pred
}

blues_me_expanded <- dplyr::bind_rows(
  expanded_blue_list
)

# Check expanded multi-environment BLUEs

print(
  blues_me_expanded %>%
    dplyr::group_by(TRAIT) %>%
    dplyr::summarise(
      N = dplyr::n(),
      N_GENOTYPES = dplyr::n_distinct(germplasmName),
      .groups = "drop"
    )
)

stopifnot(
  length(
    intersect(
      selected_2026,
      unique(blues_me_expanded$germplasmName)
    )
  ) == 0
)

saveRDS(
  blues_me_expanded,
  "results/blues_me_expanded.rds"
)

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
  "results/cv_predictions_expanded.rds"
)

write.csv(
  cv_results_expanded,
  "results/cv_results_expanded.csv",
  row.names = FALSE
)

write.csv(
  cv_summary_expanded,
  "results/cv_summary_expanded.csv",
  row.names = FALSE
)