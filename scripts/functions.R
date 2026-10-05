result_path <- function(name) file.path(config$output, name)
save_result <- function(x, name) saveRDS(x, result_path(name))
load_result <- function(name) {
  path <- result_path(name)
  if (!file.exists(path)) stop("Missing pipeline output: ", path, ". Run prepare first.")
  readRDS(path)
}
write_table <- function(x, name) {
  utils::write.csv(x, result_path(file.path("tables", name)), row.names = FALSE)
}
required_columns <- function(dat, columns, label) {
  missing <- setdiff(columns, names(dat))
  if (length(missing)) stop(label, " lacks columns: ", paste(missing, collapse = ", "))
}
clean_id <- function(x) gsub("[^A-Z0-9]", "", toupper(trimws(x)))
check_mapping <- function(dat, id) {
  x <- unique(dat[c(id, "FullSampleName")])
  x <- x[!is.na(x$FullSampleName) & nzchar(x$FullSampleName), , drop = FALSE]
  if (anyDuplicated(x[[id]])) stop("Multiple marker IDs for a genotype in ", id)
}
harmonize <- function(dat, grm) {
  ids <- rownames(grm)
  keys <- clean_id(ids)
  if (anyDuplicated(keys)) stop("Ambiguous normalized marker identifiers.")
  idx <- match(clean_id(dat$FullSampleName), keys)
  matched <- !is.na(idx)
  dat$FullSampleName[matched] <- ids[idx[matched]]
  dat
}
assert_unique <- function(dat, keys, label) {
  if (anyDuplicated(dat[keys])) stop("Duplicate records in ", label, ": ", paste(keys, collapse = ", "))
}
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(NA_real_)
  stats::cor(x[ok], y[ok])
}
model_log <- list()
fit_model <- function(..., label) {
  warnings <- character()
  arguments <- list(...)
  # sommer resolves covariance symbols in its calling frame rather than
  # consistently using the formula environment. Bind them explicitly here.
  if (!is.null(arguments$random)) for (symbol in all.vars(arguments$random)) {
    formula_env <- environment(arguments$random)
    if (!symbol %in% names(arguments$data) && exists(symbol, envir = formula_env, inherits = TRUE)) {
      assign(symbol, get(symbol, envir = formula_env, inherits = TRUE))
    }
  }
  fit <- withCallingHandlers(
    sommer::mmes(..., dateWarning = FALSE, verbose = FALSE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  convergence <- fit$convergence
  model_log[[length(model_log) + 1L]] <<- data.frame(
    model = label,
    convergence = if (is.null(convergence)) NA_character_ else paste(convergence, collapse = ";"),
    warnings = paste(unique(warnings), collapse = "; ")
  )
  if (is.logical(convergence) && any(!convergence)) stop("Model did not converge: ", label)
  fit
}

# Common fixed-genotype, random-replication phenotype model.
single_estimates <- function(dat) {
  output <- list()
  for (env in unique(dat$environment)) for (trait in config$traits) {
    x <- dat[dat$environment == env & dat$TRAIT == trait & is.finite(dat$y), ]
    if (!nrow(x)) next
    if (anyNA(x$genotype) || anyNA(x$REP)) stop("Missing genotype/replication in ", env)
    mapping <- unique(x[c("genotype", "FullSampleName")])
    check_mapping(mapping, "genotype")
    x$genotype <- factor(x$genotype)
    x$REP <- factor(x$REP)
    if (nlevels(x$REP) < 2L || nrow(x) == nlevels(x$genotype)) {
      pred <- dplyr::summarise(dplyr::group_by(x, genotype), predicted.value = mean(y), .groups = "drop")
      pred$std.error <- NA_real_
      pred$adjusted <- FALSE
    } else {
      fit <- fit_model(fixed = y ~ genotype, random = ~ REP, rcov = ~ units,
                       data = x, label = paste("phenotype", env, trait))
      pred <- as.data.frame(sommer::predict.mmes(fit, D = "genotype")$pvals)
      pred$adjusted <- TRUE
    }
    pred$genotype <- as.character(pred$genotype)
    pred <- dplyr::left_join(pred, mapping, by = "genotype")
    pred$TRAIT <- trait
    pred$environment <- env
    output[[length(output) + 1L]] <- pred[c("genotype", "FullSampleName", "TRAIT", "environment", "predicted.value", "std.error", "adjusted")]
  }
  dplyr::bind_rows(output)
}

# Preserve the existing precision-weighting approach, including mean-SE replacement.
across_estimates <- function(dat, label) {
  output <- list()
  for (trait in config$traits) {
    x <- dat[dat$TRAIT == trait & is.finite(dat$predicted.value), ]
    if (!nrow(x)) next
    assert_unique(x, c("genotype", "TRAIT", "environment"), label)
    mapping <- unique(x[c("genotype", "FullSampleName")])
    check_mapping(mapping, "genotype")
    x$y <- x$predicted.value
    x$genotype <- factor(x$genotype)
    x$environment <- factor(x$environment)
    if (nlevels(x$environment) == 1L) {
      pred <- x[c("genotype", "predicted.value", "std.error", "adjusted")]
    } else {
      args <- list(fixed = y ~ genotype, random = ~ environment, rcov = ~ units,
                   data = x, label = paste(label, trait))
      if (any(is.finite(x$std.error))) {
        se <- x$std.error
        if (any(is.finite(se) & se <= 0)) stop("Non-positive SE in ", label, " ", trait)
        se[!is.finite(se)] <- mean(se[is.finite(se)])
        args$W <- diag(1 / se^2)
      }
      fit <- do.call(fit_model, args)
      pred <- as.data.frame(sommer::predict.mmes(fit, D = "genotype")$pvals)
      pred$adjusted <- TRUE
    }
    pred$genotype <- as.character(pred$genotype)
    pred <- dplyr::left_join(pred, mapping, by = "genotype")
    pred$TRAIT <- trait
    output[[trait]] <- pred[c("genotype", "FullSampleName", "TRAIT", "predicted.value", "std.error", "adjusted")]
  }
  dplyr::bind_rows(output)
}

prediction_input <- function(dat, trait, grm) {
  x <- dat[dat$TRAIT == trait & is.finite(dat$predicted.value) &
             !is.na(dat$FullSampleName) & dat$FullSampleName %in% rownames(grm), ]
  x <- unique(x[c("FullSampleName", "predicted.value")])
  assert_unique(x, "FullSampleName", paste("prediction input", trait))
  x
}
predict_gblup <- function(train, target, grm, label) {
  assert_unique(train, "FullSampleName", paste(label, "training"))
  assert_unique(target, "FullSampleName", paste(label, "validation"))
  if (length(intersect(train$FullSampleName, target$FullSampleName))) stop("Training/validation overlap: ", label)
  if (nrow(train) < 3L || !nrow(target)) stop("Insufficient observations: ", label)
  dat <- rbind(train, target)
  dat$y <- c(train$predicted.value, rep(NA_real_, nrow(target)))
  ids <- dat$FullSampleName
  k <- grm[ids, ids, drop = FALSE]
  dat$FullSampleName <- factor(ids, levels = ids)
  fit <- fit_model(fixed = y ~ 1,
                   random = ~ sommer::vsm(sommer::ism(FullSampleName), Gu = k),
                   rcov = ~ units, data = dat, label = label)
  pred <- as.data.frame(sommer::predict.mmes(fit, D = "FullSampleName")$pvals)
  target$GBLUP <- pred$predicted.value[match(target$FullSampleName, as.character(pred$FullSampleName))]
  if (any(!is.finite(target$GBLUP))) stop("Missing/non-finite validation predictions: ", label)
  names(target)[names(target) == "predicted.value"] <- "BLUE"
  target
}

fit_heritability <- function(...) {
  notes <- character()
  fit <- withCallingHandlers(lme4::lmer(...),
    warning = function(w) { notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { notes <<- c(notes, conditionMessage(m)); invokeRestart("muffleMessage") })
  model_log[[length(model_log) + 1L]] <<- data.frame(
    model = paste("heritability", paste(deparse(fit@call$formula), collapse = " ")),
    convergence = as.character(is.null(fit@optinfo$conv$lme4$messages)),
    warnings = paste(unique(notes), collapse = "; "))
  fit
}

preparation_signature <- function() {
  list(inputs = tools::md5sum(c(config$historical, config$testing_2025, config$locations, config$markers)),
       code = tools::md5sum(c("scripts/prepare.R", "scripts/functions.R")),
       settings = config[c("traits", "maf", "heterozygosity", "unreplicated_2025")])
}
check_prepared_inputs <- function() {
  path <- result_path("preparation_manifest.rds")
  if (!file.exists(path) || !identical(readRDS(path), preparation_signature())) {
    stop("Prepared inputs are missing or stale. Run --stage prepare with the same --output first.")
  }
}
