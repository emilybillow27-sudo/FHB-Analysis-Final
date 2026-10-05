# Run from any directory: Rscript /path/to/run_analysis.R --stage all
args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg)) {
  root <- dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE)))
} else {
  root <- getwd()
}
setwd(root)
if (!identical(Sys.getenv("RENV_PROJECT"), root)) source(".Rprofile")
source("config.R")
stage <- "all"
if (length(args) %% 2L) stop("Arguments require values: --stage, --repetitions, --output")
if (length(args)) for (i in seq(1L, length(args), by = 2L)) {
  value <- args[i + 1L]
  switch(args[i],
         "--stage" = { stage <- value },
         "--repetitions" = { config$repetitions <- as.integer(value) },
         "--output" = { config$output <- value },
         stop("Unknown argument: ", args[i]))
}
stages <- c("prepare", "prediction", "cv", "summaries", "figures")
if (!stage %in% c("all", stages)) stop("Unknown stage: ", stage)
if (is.na(config$repetitions) || config$repetitions < 1L) stop("Repetitions must be a positive integer.")
packages <- c("dplyr", "tidyr", "sommer", "gaston", "ggplot2", "lme4", "e1071")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), ". Run renv::restore().")
for (dir in c(config$output, result_dir <- file.path(config$output, "tables"), file.path(config$output, "figures"))) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}
for (file in c("functions", "prepare", "prediction", "summaries", "figures")) source(file.path("scripts", paste0(file, ".R")))
source_hashes <- tools::md5sum(c("config.R", "run_analysis.R", list.files("scripts", pattern = "[.]R$", full.names = TRUE)))
started <- Sys.time()
succeeded <- FALSE
tryCatch({
  for (s in if (stage == "all") stages else stage) {
    message("Running ", s)
    if (s != "prepare") check_prepared_inputs()
    switch(s, prepare = prepare_analysis(), prediction = forward_analysis(), cv = cross_validation(),
           summaries = summary_analysis(), figures = figure_analysis())
  }
  succeeded <- TRUE
}, finally = {
  if (length(model_log)) write_table(dplyr::bind_rows(model_log), paste0("model_diagnostics_", stage, ".csv"))
  saveRDS(list(config = config, stage = stage, succeeded = succeeded, started = started, finished = Sys.time(),
               inputs = tools::md5sum(c(config$historical, config$testing_2025, config$locations, config$markers)),
               code = source_hashes,
               session = utils::sessionInfo()), result_path(paste0("run_manifest_", stage, ".rds")))
})
