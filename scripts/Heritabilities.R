# Heritabilities.R
# Author: Emily Billow
# Purpose: Estimate broad-sense heritability for the historical, 2025,
# and 2026 FHB evaluation populations.
#
# Run from the project root:
# source("scripts/Heritabilities.R")

library(dplyr)
library(tidyr)
library(lme4)
library(purrr)

traits <- c("INC", "SEV", "DON", "FDK")
traits_2026 <- traits

if (!file.exists("data/FHB_Project_Training_Data.csv")) stop("Run from the project root.")
# H2 below is an approximate entry-mean estimate using arithmetic mean
# replication/environment counts for unbalanced data, not Cullis heritability.
# Historical genotype-by-study variance can be weakly identified when most
# genotype-study cells have one observation. Retain diagnostic flags for review.
fit_notes <- character()
fit_h2_model <- function(...) {
  fit_notes <<- character()
  withCallingHandlers(lme4::lmer(...),
    warning = function(w) { fit_notes <<- c(fit_notes, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { fit_notes <<- c(fit_notes, conditionMessage(m)); invokeRestart("muffleMessage") })
}

# Locate a column using exact names or regular expressions
find_column <- function(dat, exact = character(), pattern = NULL) {
  candidates <- intersect(exact, names(dat))

  if (length(candidates) > 0) {
    return(candidates[1])
  }

  if (!is.null(pattern)) {
    matches <- grep(pattern, names(dat), ignore.case = TRUE, value = TRUE)

    if (length(matches) > 0) {
      return(matches[1])
    }
  }

  return(NA_character_)
}

# Rename a column only when a suitable source column exists
rename_trait_column <- function(dat, new_name, exact, pattern) {
  old_name <- find_column(dat, exact = exact, pattern = pattern)

  if (!is.na(old_name) && old_name != new_name) {
    names(dat)[names(dat) == old_name] <- new_name
  }

  dat
}

# Extract a variance component from an lme4 model
get_variance_component <- function(model, group_name) {
  vc <- as.data.frame(lme4::VarCorr(model))
  value <- vc$vcov[vc$grp == group_name]

  if (length(value) == 0 || is.na(value[1])) {
    return(0)
  }

  value[1]
}

# Calculate repeatability-style broad-sense heritability
calculate_h2 <- function(Vg, Ve, mean_reps, mean_environments = 1,
                         Vge = 0) {
  denominator <- Vg +
    Vge / mean_environments +
    Ve / (mean_environments * mean_reps)

  if (denominator <= 0 || is.na(denominator)) {
    return(NA_real_)
  }

  Vg / denominator
}

# Fit a single-location heritability model
calculate_single_location_h2 <- function(dat, trait, population, location) {
  df <- dat %>%
    dplyr::select(GENOTYPE, REP, y = dplyr::all_of(trait)) %>%
    dplyr::filter(is.finite(y)) %>%
    dplyr::mutate(
      GENOTYPE = as.factor(GENOTYPE),
      REP = as.factor(REP)
    )

  if (nrow(df) == 0) {
    warning("No observations for ", trait, " at ", location)
    return(NULL)
  }

  n_genotypes <- dplyr::n_distinct(df$GENOTYPE)
  n_reps <- dplyr::n_distinct(df$REP)

  if (n_genotypes < 2 || n_genotypes >= nrow(df) || n_reps < 2) {
    warning("Fewer than two genotypes for ", trait, " at ", location)
    return(NULL)
  }

  if (n_reps >= 2) {
    model <- fit_h2_model(
      y ~ 1 + (1 | GENOTYPE) + (1 | REP),
      data = df,
      REML = TRUE,
      control = lme4::lmerControl(
        optimizer = "bobyqa"
      )
    )
  } else {
    model <- fit_h2_model(
      y ~ 1 + (1 | GENOTYPE),
      data = df,
      REML = TRUE,
      control = lme4::lmerControl(
        optimizer = "bobyqa"
      )
    )
  }

  Vg <- get_variance_component(model, "GENOTYPE")
  Ve <- get_variance_component(model, "Residual")

  mean_reps <- df %>%
    dplyr::count(GENOTYPE) %>%
    dplyr::summarise(mean_reps = mean(n)) %>%
    dplyr::pull(mean_reps)

  H2 <- if (n_reps >= 2) {
    calculate_h2(
      Vg = Vg,
      Ve = Ve,
      mean_reps = mean_reps
    )
  } else {
    NA_real_
  }

  data.frame(
    Trait = trait,
    Population = population,
    Location = location,
    n_observations = nrow(df),
    n_genotypes = n_genotypes,
    mean_reps_per_genotype = mean_reps,
    Vg = Vg,
    Vge = NA_real_,
    Ve = Ve,
    H2 = H2,
    singular_fit = lme4::isSingular(model, tol = 1e-4),
    convergence_warning = paste(model@optinfo$conv$lme4$messages, collapse = "; "),
    model_notes = paste(unique(fit_notes), collapse = "; "),
    estimator = "Approximate entry-mean H2"
  )
}

# Read the original 2025 testing data
testing_2025 <- read.csv(
  "data/FHB_Project_Testing_Data_2025.csv",
  check.names = FALSE, fileEncoding = "UTF-8-BOM"
)

# Standardize 2025 identifiers and traits
testing_2025 <- testing_2025 %>%
  dplyr::mutate(
    GENOTYPE = as.character(ID),
    REP = as.factor(REP)
  )

for (trait in traits) {
  if (!trait %in% names(testing_2025)) {
    testing_2025 <- rename_trait_column(
      testing_2025,
      trait,
      exact = trait,
      pattern = trait
    )
  }
}

# Read the original historical training data
training <- read.csv(
  "data/FHB_Project_Training_Data.csv",
  check.names = FALSE, fileEncoding = "UTF-8-BOM"
)

# Rename historical traits
training <- training %>%
  rename_trait_column(
    "DON",
    exact = "FHB.DON.content...ppm.CO_321.0001154",
    pattern = "DON.content"
  ) %>%
  rename_trait_column(
    "INC",
    exact = "FHB.incidence.....CO_321.0001149",
    pattern = "^FHB[.]incidence"
  ) %>%
  rename_trait_column(
    "SEV",
    exact = "FHB.severity.....CO_321.0001440",
    pattern = "severity"
  ) %>%
  rename_trait_column(
    "FDK",
    exact = "FHB.grain.incidence.....CO_321.0001155",
    pattern = "grain.incidence|FDK"
  )

# Confirm historical identifiers
if (!"germplasmName" %in% names(training)) {
  stop("The historical training data require a germplasmName column.")
}

if (!"studyName" %in% names(training)) {
  stop("The historical training data require a studyName column.")
}

if (!"replicate" %in% names(training)) {
  stop("The historical training data require a replicate column.")
}

training <- training %>%
  dplyr::mutate(
    GENOTYPE = as.character(germplasmName),
    STUDY = as.character(studyName),
    REP = as.factor(replicate)
  )

# Analyze replicated nursery components separately, as in the main pipeline.
h2_2025 <- purrr::map_dfr(setdiff(unique(testing_2025$SUB_NURNAME), c("DH", "Topcross", "Q Qual AYN")), function(nursery) {
  dat <- testing_2025[testing_2025$SUB_NURNAME == nursery, ]
  purrr::map_dfr(intersect(traits, names(dat)), function(trait) {
    calculate_single_location_h2(dat, trait, paste("2025 testing", nursery), "UIUC")
  })
})

# Calculate historical across-environment heritability
calculate_h2_historical <- function(dat, trait) {
  df <- dat %>%
    dplyr::select(GENOTYPE, STUDY, REP, y = dplyr::all_of(trait)) %>%
    dplyr::filter(is.finite(y)) %>%
    dplyr::mutate(
      GENOTYPE = as.factor(GENOTYPE),
      STUDY = as.factor(STUDY),
      REP = as.factor(REP)
    )

  if (nrow(df) == 0) {
    warning("No historical observations for ", trait)
    return(NULL)
  }

  model <- fit_h2_model(
    y ~ 1 +
      (1 | GENOTYPE) +
      (1 | STUDY) +
      (1 | GENOTYPE:STUDY) +
      (1 | STUDY:REP),
    data = df,
    REML = TRUE,
    control = lme4::lmerControl(
      optimizer = "bobyqa"
    )
  )

  Vg <- get_variance_component(model, "GENOTYPE")
  Vge <- get_variance_component(model, "GENOTYPE:STUDY")
  Ve <- get_variance_component(model, "Residual")

  genotype_environment_counts <- df %>%
    dplyr::group_by(GENOTYPE, STUDY) %>%
    dplyr::summarise(n_obs = dplyr::n(), .groups = "drop")

  mean_environments <- genotype_environment_counts %>%
    dplyr::count(GENOTYPE) %>%
    dplyr::summarise(value = mean(n)) %>%
    dplyr::pull(value)

  mean_reps <- mean(genotype_environment_counts$n_obs)

  H2 <- calculate_h2(
    Vg = Vg,
    Vge = Vge,
    Ve = Ve,
    mean_reps = mean_reps,
    mean_environments = mean_environments
  )

  data.frame(
    Trait = trait,
    Population = "Historical training",
    Location = "Multiple environments",
    n_observations = nrow(df),
    n_genotypes = dplyr::n_distinct(df$GENOTYPE),
    mean_reps_per_genotype = mean_reps,
    Vg = Vg,
    Vge = Vge,
    Ve = Ve,
    H2 = H2,
    singular_fit = lme4::isSingular(model, tol = 1e-4),
    convergence_warning = paste(model@optinfo$conv$lme4$messages, collapse = "; "),
    model_notes = paste(unique(fit_notes), collapse = "; "),
    estimator = "Approximate entry-mean H2"
  )
}

historical_traits <- intersect(traits, names(training))

h2_historical <- purrr::map_dfr(
  historical_traits,
  ~ calculate_h2_historical(training, .x)
)

# Read the original 2026 datasets
standardize_2026 <- function(file_name, location) {
  dat <- read.csv(
    file_name,
    check.names = FALSE
  )

  required_columns <- c("ID", "REP", "Incidence", "Severity")
  missing_columns <- setdiff(required_columns, names(dat))

  if (length(missing_columns) > 0) {
    stop(
      location,
      " is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  dat <- dat %>%
    dplyr::rename(
      INC = Incidence,
      SEV = Severity
    ) %>%
    dplyr::mutate(
      GENOTYPE = as.character(ID),
      REP = as.factor(REP)
    ) %>% dplyr::filter(GENOTYPE != "Triticale")

  dat %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(intersect(traits_2026, names(.))),
        as.numeric
      )
    )
}

# Use all available wheat entries at each location; absent traits are skipped.
h2_2026 <- purrr::map_dfr(c("UNL", "UIUC", "SDSU"), function(location) {
  dat <- standardize_2026(paste0("data/FHB_Project_Testing_Data_", location, "_2026.csv"), location)
  purrr::map_dfr(intersect(traits, names(dat)), function(trait) {
    calculate_single_location_h2(dat, trait, "2026 testing", location)
  })
})

# Combine all results
heritability_results <- dplyr::bind_rows(
  h2_historical,
  h2_2025,
  h2_2026
)

# Create the output directory
dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

# Save the results
write.csv(
  heritability_results,
  "results/tables/heritabilities.csv",
  row.names = FALSE
)

# Print the results
print(heritability_results)

# Print model warning summary
cat("\nSingular fits are flagged in the singular_fit column.\n")
write.csv(heritability_results[, c("Trait", "Population", "Location", "singular_fit", "convergence_warning", "model_notes")],
          "results/tables/heritability_diagnostics.csv", row.names = FALSE)
