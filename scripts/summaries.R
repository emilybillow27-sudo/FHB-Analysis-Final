# Shared descriptive summaries and existing heritability model specifications.
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
  df <- dat |>
    dplyr::select(GENOTYPE, REP, y = dplyr::all_of(trait)) |>
    dplyr::filter(!is.na(y)) |>
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

  if (n_genotypes < 2 || n_genotypes >= nrow(df)) {
    warning("Fewer than two genotypes for ", trait, " at ", location)
    return(NULL)
  }

  if (n_reps >= 2) {
    model <- fit_heritability(
      y ~ 1 + (1 | GENOTYPE) + (1 | REP),
      data = df,
      REML = TRUE,
      control = lme4::lmerControl(
        optimizer = "bobyqa"
      )
    )
  } else {
    model <- fit_heritability(
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

  mean_reps <- df |>
    dplyr::count(GENOTYPE) |>
    dplyr::summarise(mean_reps = mean(n)) |>
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
    convergence_warning = paste(model@optinfo$conv$lme4$messages, collapse = "; ")
  )
}

calculate_h2_historical <- function(dat, trait) {
  df <- dat |>
    dplyr::select(GENOTYPE, STUDY, REP, y = dplyr::all_of(trait)) |>
    dplyr::filter(!is.na(y)) |>
    dplyr::mutate(
      GENOTYPE = as.factor(GENOTYPE),
      STUDY = as.factor(STUDY),
      REP = as.factor(REP)
    )

  if (nrow(df) == 0) {
    warning("No historical observations for ", trait)
    return(NULL)
  }

  model <- fit_heritability(
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

  genotype_environment_counts <- df |>
    dplyr::group_by(GENOTYPE, STUDY) |>
    dplyr::summarise(n_obs = dplyr::n(), .groups = "drop")

  mean_environments <- genotype_environment_counts |>
    dplyr::count(GENOTYPE) |>
    dplyr::summarise(value = mean(n)) |>
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
    convergence_warning = paste(model@optinfo$conv$lme4$messages, collapse = "; ")
  )
}


summary_analysis <- function() {
  raw <- load_result("raw.rds")
  historical <- raw$historical
  historical$GENOTYPE <- historical$genotype
  historical$STUDY <- historical$environment
  historical$REP <- factor(historical$REP)
  populations <- c(list(Historical = raw$historical, `2025 UIUC` = raw$testing), raw$locations)
  summaries <- list()
  for (name in names(populations)) {
    x <- long_phenotypes(populations[[name]])
    summary <- dplyr::summarise(dplyr::group_by(x, TRAIT),
      n = sum(is.finite(y)), mean = mean(y[is.finite(y)]), sd = stats::sd(y[is.finite(y)]),
      minimum = if (any(is.finite(y))) min(y[is.finite(y)]) else NA_real_,
      maximum = if (any(is.finite(y))) max(y[is.finite(y)]) else NA_real_,
      skewness = if (sum(is.finite(y)) >= 3L && stats::sd(y[is.finite(y)]) > 0)
        e1071::skewness(y[is.finite(y)], type = 2) else NA_real_, .groups = "drop")
    summary$Population <- name
    summaries[[name]] <- summary
  }
  write_table(dplyr::bind_rows(summaries), "phenotype_summary.csv")
  heritability <- lapply(intersect(config$traits, names(historical)), function(trait) calculate_h2_historical(historical, trait))
  for (name in setdiff(names(populations), "Historical")) {
    x <- populations[[name]]
    x$GENOTYPE <- x$genotype
    for (trait in intersect(config$traits, names(x))) {
      if (!any(is.finite(x[[trait]]))) next
      heritability[[length(heritability) + 1L]] <- calculate_single_location_h2(
        x, trait, if (name == "2025 UIUC") "2025 testing" else "2026 testing",
        if (name == "2025 UIUC") "UIUC" else name)
    }
  }
  write_table(dplyr::bind_rows(heritability), "heritabilities.csv")
}
