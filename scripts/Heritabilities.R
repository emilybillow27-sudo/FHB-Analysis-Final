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

bootstrap_n <- as.integer(Sys.getenv("FHB_H2_BOOTSTRAPS", "500"))
if (is.na(bootstrap_n) || bootstrap_n < 100L) stop("Use at least 100 bootstrap simulations.")
model_fits <- list()

traits <- c("INC", "SEV", "DON", "FDK")
traits_2026 <- traits

if (!file.exists("data/FHB_Project_Training_Data.csv")) stop("Run from the project root.")
# H2 below is an approximate entry-mean estimate using arithmetic mean
# replication/environment counts for unbalanced data, not Cullis heritability.
# Historical genotype-by-study variance can be weakly identified when most
# genotype-study cells have one observation. Retain diagnostic flags for review.
fit_notes <- character()
fit_h2_model <- function(formula, data, REML = TRUE, control = lme4::lmerControl(optimizer = "bobyqa")) {
  fit_notes <<- character()
  model <- withCallingHandlers(lme4::lmer(formula = formula, data = data, REML = REML, control = control),
    warning = function(w) { fit_notes <<- c(fit_notes, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { fit_notes <<- c(fit_notes, conditionMessage(m)); invokeRestart("muffleMessage") })
  model@call$formula <- formula
  model@call$data <- data
  model@call$REML <- REML
  model@call$control <- quote(lme4::lmerControl(optimizer = "bobyqa"))
  model

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

# Report every fitted variance component and the source of boundary fits.
model_diagnostics <- function(model, df, historical = FALSE) {
  vc <- as.data.frame(lme4::VarCorr(model))
  random <- vc[vc$grp != "Residual", ]
  near_zero <- random$grp[random$sdcor <= 1e-4 * sigma(model)]
  cells <- if (historical) interaction(df$GENOTYPE, df$STUDY, drop = TRUE) else df$GENOTYPE
  counts <- table(cells)
  data.frame(
    Vrep = get_variance_component(model, if (historical) "STUDY:REP" else "REP"),
    Vstudy = if (historical) get_variance_component(model, "STUDY") else NA_real_,
    boundary_components = paste(near_zero, collapse = "; "),
    singleton_cells = sum(counts == 1),
    replicated_cells = sum(counts > 1),
    mean_environments = if (historical) length(counts) / nlevels(df$GENOTYPE) else 1,
    review_status = if ("GENOTYPE" %in% near_zero) "Genotype variance at boundary" else
      if (length(near_zero)) "Other random-effect variance at boundary" else
      if (length(fit_notes)) "Review model warnings" else "No fitting flags"
  )
}

# Generalized heritability uses the full prediction-error covariance, including
# uncertainty in fixed effects. Undefined when genotype variance is zero.
cullis_h2 <- function(model) {
  vg <- get_variance_component(model, 'GENOTYPE')
  if(vg <= 1e-8 * sigma(model)^2) return(NA_real_)
  X <- lme4::getME(model,'X'); A <- lme4::getME(model,'Z') %*% lme4::getME(model,'Lambda')
  W <- cbind(Matrix::Matrix(X,sparse=TRUE),A); p <- ncol(X)
  Q <- crossprod(W)+Matrix::Diagonal(ncol(W),c(rep(0,p),rep(1,ncol(A))))
  groups <- names(lme4::getME(model,'cnms')); k <- match('GENOTYPE',groups); gp <- lme4::getME(model,'Gp')
  ids <- p + seq.int(gp[k]+1,gp[k+1]); n <- length(ids)
  rhs <- Matrix::sparseMatrix(i=ids,j=seq_len(n),x=1,dims=c(ncol(W),n))
  C <- as.matrix(Matrix::solve(Q,rhs)[ids,,drop=FALSE])
  1-(sum(diag(C))-sum(C)/n)/(n-1)
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

  model_fits[[paste(population, location, trait, sep = "|")]] <<- model

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

  cbind(model_diagnostics(model, df), data.frame(
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
  ))
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

  cbind(model_diagnostics(model, df, historical = TRUE), data.frame(
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
  ))
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

# Sensitivity analysis for the replicated validation population, excluding checks
# and single-replicate entries. Keep the full-nursery estimates above unchanged.
validation_ids <- read.csv("data/FHB_Project_Testing_Data_UNL_2026.csv", check.names = FALSE)
validation_ids <- unique(validation_ids$ID[!is.na(validation_ids$CATEGORY) & validation_ids$CATEGORY == "Duplicate"])
if (length(validation_ids) != 60L) stop("Expected 60 validation genotypes.")
h2_validation <- purrr::map_dfr(c("UNL", "UIUC", "SDSU"), function(location) {
  dat <- standardize_2026(paste0("data/FHB_Project_Testing_Data_", location, "_2026.csv"), location)
  dat <- dat[dat$GENOTYPE %in% validation_ids, ]
  purrr::map_dfr(intersect(traits, names(dat)), function(trait) {
    calculate_single_location_h2(dat, trait, "2026 validation lines", location)
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

# Pipeline-stage table folders
dir.create("results/tables/heritability", recursive=TRUE, showWarnings=FALSE)
write.csv(
  heritability_results,
  "results/tables/heritability/heritabilities.csv",
  row.names = FALSE
)

write.csv(h2_validation, "results/tables/heritability/heritability_validation_sensitivity.csv", row.names = FALSE)

# Print the results
print(heritability_results)

# Print model warning summary
cat("\nSingular fits are flagged in the singular_fit column.\n")
write.csv(heritability_results[, c("Trait", "Population", "Location", "Vg", "Vge", "Vrep", "Vstudy", "Ve", "boundary_components", "singleton_cells", "replicated_cells", "mean_environments", "review_status", "singular_fit", "convergence_warning", "model_notes")],
          "results/tables/heritability/heritability_diagnostics.csv", row.names = FALSE)

# Primary reporting population: replicated 2025 nurseries and 2026 validation lines.
# Preserve the full-nursery and historical exploratory outputs above.
heritability_primary <- dplyr::bind_rows(h2_2025, h2_validation)
set.seed(20261005)
for (i in seq_len(nrow(heritability_primary))) {
  row <- heritability_primary[i, ]
  model <- model_fits[[paste(row$Population, row$Location, row$Trait, sep = "|")]]
  reps <- row$mean_reps_per_genotype
  message("Bootstrapping ", row$Population, " ", row$Location, " ", row$Trait)
  statistic <- function(fit) {
    vg <- get_variance_component(fit, "GENOTYPE")
    ve <- get_variance_component(fit, "Residual")
    bad <- fit@optinfo$conv$lme4$messages
    bad <- bad[!grepl("boundary.*singular", bad)]
    c(H2 = calculate_h2(vg, ve, reps), singular = as.numeric(lme4::isSingular(fit)),
      convergence_flag = as.numeric(length(bad) > 0))
  }
  boot <- lme4::bootMer(model, FUN = statistic, nsim = bootstrap_n,
                        use.u = FALSE, type = "parametric")
  valid <- is.finite(boot$t[, 1]) & boot$t[, 3] == 0
  limits <- if (sum(valid) >= 0.9 * bootstrap_n) quantile(boot$t[valid, 1], c(.025, .975)) else c(NA, NA)
  heritability_primary$H2_lower[i] <- limits[1]
  heritability_primary$H2_upper[i] <- limits[2]
  heritability_primary$bootstrap_requested[i] <- bootstrap_n
  heritability_primary$bootstrap_valid[i] <- sum(valid)
  heritability_primary$bootstrap_singular_fraction[i] <- mean(boot$t[valid, 2])
  heritability_primary$uncertainty_method[i] <- if ("GENOTYPE" %in% strsplit(row$boundary_components, "; ")[[1]])
    "Parametric bootstrap percentile range; boundary estimate, nominal coverage uncertain" else
    "Approximate 95% parametric bootstrap percentile interval"
}
primary_columns <- c("Trait", "Population", "Location", "n_observations", "n_genotypes",
                     "H2", "H2_lower", "H2_upper", "boundary_components",
                     "uncertainty_method", "bootstrap_requested", "bootstrap_valid",
                     "bootstrap_singular_fraction", "estimator")
write.csv(heritability_primary[, primary_columns], "results/tables/heritability/heritability_primary.csv", row.names = FALSE)

# Historical design audit and exploratory block sensitivity.
# Do not select a model merely to remove a warning. Confirm these block labels
# against the historical field layouts before using historical H2 in the thesis.
historical_design <- list()
historical_sensitivity <- list()
for (trait in historical_traits) {
  df <- training[is.finite(training[[trait]]), ]
  df$y <- df[[trait]]
  cells <- dplyr::count(df, GENOTYPE, STUDY)
  coverage <- dplyr::count(cells, GENOTYPE)
  historical_design[[trait]] <- data.frame(
    Trait = trait, n_observations = nrow(df), n_genotypes = nrow(coverage),
    n_studies = dplyr::n_distinct(df$STUDY),
    genotypes_one_study = sum(coverage$n == 1),
    genotypes_multiple_studies = sum(coverage$n > 1),
    singleton_cells = sum(cells$n == 1), replicated_cells = sum(cells$n > 1),
    unique_plot_units = dplyr::n_distinct(df$observationUnitDbId))
  df$BLOCK <- ifelse(df$studyDesign == "RCBD" & !is.na(df$blockNumber),
                     as.character(df$blockNumber), as.character(df$REP))
  for (specification in c("Recorded replicates", "Recorded RCBD blocks")) {
    formula <- if (specification == "Recorded replicates")
      y ~ 1 + (1 | GENOTYPE) + (1 | STUDY) + (1 | GENOTYPE:STUDY) + (1 | STUDY:REP) else
      y ~ 1 + (1 | GENOTYPE) + (1 | STUDY) + (1 | GENOTYPE:STUDY) + (1 | STUDY:BLOCK)
    model <- fit_h2_model(formula, data = df)
    historical_sensitivity[[length(historical_sensitivity) + 1L]] <- data.frame(
      Trait = trait, block_specification = specification,
      Vg = get_variance_component(model, "GENOTYPE"),
      Vge = get_variance_component(model, "GENOTYPE:STUDY"), Ve = sigma(model)^2,
      H2_Cullis = cullis_h2(model), model_notes = paste(fit_notes, collapse = "; "),
      status = "Exploratory: verify historical blocking before thesis use")
  }
}
write.csv(dplyr::bind_rows(historical_design), "results/tables/heritability/historical_heritability_design.csv", row.names = FALSE)
write.csv(dplyr::bind_rows(historical_sensitivity), "results/tables/heritability/heritability_historical_sensitivity.csv", row.names = FALSE)
