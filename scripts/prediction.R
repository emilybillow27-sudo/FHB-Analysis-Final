forward_analysis <- function() {
  grm <- load_result("GRM.rds")
  historical <- load_result("historical_me.rds")
  expanded <- load_result("expanded_me.rds")
  selected <- load_result("validation_ids.rds")
  testing <- load_result("testing_se.rds")
  testing <- dplyr::summarise(dplyr::group_by(testing, FullSampleName, TRAIT),
                             predicted.value = mean(predicted.value), .groups = "drop")
  evaluation <- load_result("evaluation_me.rds")
  by_location <- load_result("evaluation_se.rds")
  output <- list()
  for (trait in config$traits) {
    target <- prediction_input(testing, trait, grm)
    if (nrow(target)) {
      train <- prediction_input(historical, trait, grm)
      train <- train[!train$FullSampleName %in% target$FullSampleName, ]
      pred <- predict_gblup(train, target, grm, paste("2025", trait))
      pred$YEAR <- 2025L
      pred$environment <- "UIUC"
      pred$TRAINING_POPULATION <- "Historical"
      pred$TRAIT <- trait
      output[[length(output) + 1L]] <- pred
    }
    for (loc in c("Across locations", unique(by_location$environment))) {
      observed <- if (loc == "Across locations") evaluation else by_location[by_location$environment == loc, ]
      observed <- observed[observed$FullSampleName %in% selected, ]
      target <- prediction_input(observed, trait, grm)
      if (!nrow(target)) next
      for (scenario in c("Historical", "Expanded")) {
        training <- if (scenario == "Historical") historical else expanded
        train <- prediction_input(training, trait, grm)
        # Exclude every validation genotype, even when its phenotype is unavailable for this trait.
        train <- train[!train$FullSampleName %in% selected, ]
        pred <- predict_gblup(train, target, grm, paste("2026", loc, scenario, trait))
        pred$YEAR <- 2026L
        pred$environment <- loc
        pred$TRAINING_POPULATION <- scenario
        pred$TRAIT <- trait
        output[[length(output) + 1L]] <- pred
      }
    }
  }
  predictions <- dplyr::bind_rows(output)
  # Both 2026 scenarios must contain the same observations and usable predictions.
  pairs <- predictions[predictions$YEAR == 2026L, ]
  a <- pairs[pairs$TRAINING_POPULATION == "Historical", c("FullSampleName", "environment", "TRAIT", "BLUE")]
  b <- pairs[pairs$TRAINING_POPULATION == "Expanded", c("FullSampleName", "environment", "TRAIT", "BLUE")]
  if (!isTRUE(all.equal(dplyr::arrange(a, environment, TRAIT, FullSampleName),
                       dplyr::arrange(b, environment, TRAIT, FullSampleName), check.attributes = FALSE))) {
    stop("Historical/expanded validation observations differ.")
  }
  summary <- dplyr::summarise(dplyr::group_by(predictions, YEAR, environment, TRAINING_POPULATION, TRAIT),
                             N = dplyr::n(), PREDICTION_ABILITY = safe_cor(BLUE, GBLUP), .groups = "drop")
  comparison <- tidyr::pivot_wider(summary[summary$YEAR == 2026L, ],
                                   names_from = TRAINING_POPULATION, values_from = c(N, PREDICTION_ABILITY))
  comparison$CHANGE_IN_PREDICTION_ABILITY <- comparison$PREDICTION_ABILITY_Expanded - comparison$PREDICTION_ABILITY_Historical
  save_result(predictions, "forward_predictions.rds")
  write_table(predictions, "forward_predictions.csv")
  write_table(summary, "forward_prediction_summary.csv")
  write_table(comparison, "forward_prediction_comparison_2026.csv")
}

cross_validation <- function() {
  grm <- load_result("GRM.rds")
  output <- list()
  set.seed(config$seed)
  for (scenario in c("Historical", "Expanded")) {
    phenotypes <- load_result(if (scenario == "Historical") "historical_me.rds" else "expanded_me.rds")
    for (trait in config$traits) {
      x <- prediction_input(phenotypes, trait, grm)
      if (nrow(x) < config$folds + 3L) stop("Insufficient CV population: ", scenario, " ", trait)
      for (repetition in seq_len(config$repetitions)) {
        message("CV ", scenario, " ", trait, " ", repetition, "/", config$repetitions)
        folds <- sample(rep(seq_len(config$folds), length.out = nrow(x)))
        for (fold in seq_len(config$folds)) {
          pred <- predict_gblup(x[folds != fold, ], x[folds == fold, ], grm,
                                paste("CV", scenario, trait, repetition, fold))
          pred$TRAIT <- trait
          pred$TRAINING_POPULATION <- scenario
          pred$REPETITION <- repetition
          pred$FOLD <- fold
          output[[length(output) + 1L]] <- pred
        }
      }
    }
  }
  predictions <- dplyr::bind_rows(output)
  assert_unique(predictions, c("TRAINING_POPULATION", "TRAIT", "REPETITION", "FullSampleName"), "CV predictions")
  results <- dplyr::summarise(dplyr::group_by(predictions, TRAINING_POPULATION, TRAIT, REPETITION),
                             N = dplyr::n(), PREDICTION_ABILITY = safe_cor(BLUE, GBLUP), .groups = "drop")
  summary <- dplyr::summarise(dplyr::group_by(results, TRAINING_POPULATION, TRAIT),
                             N_REPETITIONS = dplyr::n(),
                             N_VALID_CORRELATIONS = sum(is.finite(PREDICTION_ABILITY)),
                             MEAN_PREDICTION_ABILITY = mean(PREDICTION_ABILITY, na.rm = TRUE),
                             SD_PREDICTION_ABILITY = stats::sd(PREDICTION_ABILITY, na.rm = TRUE), .groups = "drop")
  save_result(predictions, "cv_predictions.rds")
  write_table(results, "cv_results.csv")
  write_table(summary, "cv_summary.csv")
}
