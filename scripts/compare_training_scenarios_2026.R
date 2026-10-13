# Paired comparison of 2026 genomic prediction ability.
# Run from the FHB Analysis Final project folder:
# source("scripts/compare_training_scenarios_2026.R")
# Base R only. Existing prediction and analysis files are not modified.

prediction_file <- "results/tables/forward_prediction/2026/forward_predictions_2026.csv"
output_dir <- "results/tables/training_scenario_comparison_2026"
table_dir <- "results/tables/training_scenario_comparison_2026"
figure_dir <- "results/figures/training_scenario_comparison_2026"
n_boot <- 10000L
seed <- 20261009L
alpha <- 0.05
trait_order <- c("DON", "FDK", "INC", "SEV")

# Optional CSV with columns FullSampleName and Family (one row per genotype).
# Set this to a file path to resample whole families rather than genotypes.
# Family assignments must come from your breeding records, not ID prefixes.
family_file <- NULL

pred <- read.csv(prediction_file, stringsAsFactors = FALSE)
required <- c("FullSampleName", "BLUE", "GEBV", "TRAIT", "TRAINING_POPULATION")
if (!all(required %in% names(pred))) stop("Prediction file is missing required columns.")
if (!is.numeric(pred$BLUE) || !is.numeric(pred$GEBV)) stop("BLUE and GEBV must be numeric.")
pred <- pred[pred$TRAINING_POPULATION %in% c("Historical", "Expanded"), ]
if (anyNA(pred[c("FullSampleName", "TRAIT")]) || any(pred$FullSampleName == "")) {
  stop("Missing genotype or trait identifiers.")
}
if (anyDuplicated(pred[c("FullSampleName", "TRAIT", "TRAINING_POPULATION")])) {
  stop("Duplicate genotype-trait-scenario rows; resolve these before comparison.")
}

hist <- pred[pred$TRAINING_POPULATION == "Historical",
             c("FullSampleName", "TRAIT", "BLUE", "GEBV")]
expanded <- pred[pred$TRAINING_POPULATION == "Expanded",
                 c("FullSampleName", "TRAIT", "BLUE", "GEBV")]
paired <- merge(hist, expanded, by = c("FullSampleName", "TRAIT"),
                suffixes = c("_Historical", "_Expanded"), all = TRUE)
both_blue <- is.finite(paired$BLUE_Historical) & is.finite(paired$BLUE_Expanded)
if (any(abs(paired$BLUE_Historical[both_blue] - paired$BLUE_Expanded[both_blue]) >
        1e-8 * pmax(1, abs(paired$BLUE_Historical[both_blue])))) {
  stop("Observed BLUEs differ between scenarios for the same genotype and trait.")
}
paired$BLUE <- paired$BLUE_Historical

clustered <- !is.null(family_file)
if (clustered) {
  families <- read.csv(family_file, stringsAsFactors = FALSE)
  if (!all(c("FullSampleName", "Family") %in% names(families))) {
    stop("Family file must contain FullSampleName and Family.")
  }
  if (anyDuplicated(families$FullSampleName)) stop("Duplicate IDs in family file.")
  paired$Family <- families$Family[match(paired$FullSampleName, families$FullSampleName)]
}

safe_cor <- function(x, y) {
  if (length(x) < 4L || !is.finite(sd(x)) || !is.finite(sd(y)) ||
      sd(x) == 0 || sd(y) == 0) return(NA_real_)
  cor(x, y, method = "pearson")
}

set.seed(seed)
summary_list <- draws_list <- used_list <- list()
for (trait in trait_order) {
  all_dat <- paired[paired$TRAIT == trait, ]
  complete <- with(all_dat, is.finite(BLUE_Historical) & is.finite(BLUE_Expanded) &
                     is.finite(GEBV_Historical) & is.finite(GEBV_Expanded))
  dat <- all_dat[complete, ]
  n <- nrow(dat)
  if (n < 4L) stop("Fewer than four complete pairs for ", trait)
  if (clustered && (anyNA(dat$Family) || any(trimws(dat$Family) == ""))) {
    stop("Missing family assignments for complete pairs in ", trait)
  }
  groups <- if (clustered) split(seq_len(n), dat$Family) else as.list(seq_len(n))
  n_units <- length(groups)
  if (n_units < 4L) stop("Fewer than four resampling units for ", trait)
  if (clustered && n_units < 20L) warning("Few families for ", trait, "; intervals may be unstable.")
  r_hist <- safe_cor(dat$BLUE, dat$GEBV_Historical)
  r_exp <- safe_cor(dat$BLUE, dat$GEBV_Expanded)
  r_predictions <- safe_cor(dat$GEBV_Historical, dat$GEBV_Expanded)
  if (any(!is.finite(c(r_hist, r_exp)))) stop("Undefined observed correlation for ", trait)

  # Duplicate selected families contribute all their rows each time selected.
  draws <- replicate(n_boot, {
    idx <- unlist(groups[sample.int(n_units, n_units, replace = TRUE)], use.names = FALSE)
    rh <- safe_cor(dat$BLUE[idx], dat$GEBV_Historical[idx])
    re <- safe_cor(dat$BLUE[idx], dat$GEBV_Expanded[idx])
    c(Historical = rh, Expanded = re, Difference = re - rh)
  })
  valid <- apply(draws, 2, function(x) all(is.finite(x)))
  if (sum(valid) < 0.99 * n_boot) {
    stop("More than 1% undefined bootstrap samples for ", trait, "; inspect data/cluster structure.")
  }
  if (any(!valid)) warning(sum(!valid), " undefined bootstrap samples omitted for ", trait)
  ci <- function(row, a) unname(quantile(draws[row, valid], c(a / 2, 1 - a / 2)))
  diff_ci <- ci("Difference", alpha)
  # Approximate familywise coverage across all four traits via Bonferroni CIs.
  familywise_ci <- ci("Difference", alpha / length(trait_order))
  hist_ci <- ci("Historical", alpha)
  exp_ci <- ci("Expanded", alpha)
  summary_list[[trait]] <- data.frame(
    TRAIT = trait, N_PAIRED = n, N_EXCLUDED_OR_UNPAIRED = nrow(all_dat) - n,
    RESAMPLING_UNIT = if (clustered) "Family" else "Genotype", N_UNITS = n_units,
    R_HISTORICAL = r_hist, HISTORICAL_CI95_LOW = hist_ci[1], HISTORICAL_CI95_HIGH = hist_ci[2],
    R_EXPANDED = r_exp, EXPANDED_CI95_LOW = exp_ci[1], EXPANDED_CI95_HIGH = exp_ci[2],
    DELTA_R = r_exp - r_hist, DELTA_CI95_LOW = diff_ci[1], DELTA_CI95_HIGH = diff_ci[2],
    DELTA_CI98_75_LOW = familywise_ci[1], DELTA_CI98_75_HIGH = familywise_ci[2],
    CI95_EXCLUDES_ZERO = diff_ci[1] > 0 || diff_ci[2] < 0,
    FAMILYWISE_CI_EXCLUDES_ZERO = familywise_ci[1] > 0 || familywise_ci[2] < 0,
    R_BETWEEN_PREDICTIONS = r_predictions,
    BOOTSTRAP_VALID = sum(valid)
  )
  draws_list[[trait]] <- data.frame(TRAIT = trait, ITERATION = seq_len(n_boot),
                                   t(draws), VALID = valid, row.names = NULL)
  used_list[[trait]] <- dat
  message("Completed ", trait, ": ", n, " paired genotypes.")
}
comparison <- do.call(rbind, summary_list)
rownames(comparison) <- NULL
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(comparison, file.path(output_dir, "comparison_summary.csv"), row.names = FALSE, na = "")
write.csv(do.call(rbind, used_list), file.path(output_dir, "paired_genotype_predictions.csv"),
          row.names = FALSE, na = "")
write.csv(do.call(rbind, draws_list), file.path(output_dir, "bootstrap_draws.csv"),
          row.names = FALSE, na = "")

# Thesis table and figure: bootstrap only, with no text-file exports.
for (directory in c(table_dir, figure_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}
trait_labels <- c(DON = "DON", FDK = "FDK", INC = "Incidence", SEV = "Severity")
fmt <- function(x) sprintf("%.3f", x)
ci_label <- function(low, high) paste0("[", fmt(low), ", ", fmt(high), "]")
thesis_table <- data.frame(
  Trait = unname(trait_labels[comparison$TRAIT]),
  N = comparison$N_PAIRED,
  `Historical r` = fmt(comparison$R_HISTORICAL),
  `Expanded r` = fmt(comparison$R_EXPANDED),
  `Difference in r` = fmt(comparison$DELTA_R),
  `Adjusted bootstrap CI` = ci_label(comparison$DELTA_CI98_75_LOW, comparison$DELTA_CI98_75_HIGH),
  check.names = FALSE
)
unit_label <- if (clustered) "families" else "genotypes"
caption <- paste0(
  "Prediction ability is the Pearson correlation between GEBVs and observed BLUEs. ",
  "Difference = expanded minus historical. Intervals are ",
  format(100 * (1 - alpha / length(trait_order)), trim = TRUE),
  "% paired percentile bootstrap intervals (Bonferroni adjustment for four traits; ",
  "approximate 95% familywise coverage), based on ", format(n_boot, big.mark = ","),
  " resamples of ", unit_label, ". All adjusted intervals include zero. ",
  "Uncertainty is conditional on fitted models and estimated BLUEs."
)
write.csv(thesis_table, file.path(table_dir, "training_scenario_comparison_2026.csv"),
          row.names = FALSE)
# Table PNG using the same values as the formatted CSV.
draw_table <- function() {
  par(mar = rep(0, 4), family = "serif")
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0.04, 0.91, "Comparison of training populations in 2026", adj = 0, cex = 1.2)
  xpos <- c(0.04, 0.24, 0.37, 0.51, 0.66, 0.88)
  headers <- c("Trait", "N", "Historical r", "Expanded r", "Difference", "Adjusted CI")
  segments(0.04, 0.81, 0.97, 0.81, lwd = 1.2)
  for (j in seq_along(xpos)) text(xpos[j], 0.76, headers[j], adj = if (j == 1) 0 else 0.5, cex = 0.85)
  segments(0.04, 0.71, 0.97, 0.71)
  for (i in seq_len(nrow(thesis_table))) {
    for (j in seq_along(xpos)) text(xpos[j], 0.66 - (i-1)*0.085,
      as.character(thesis_table[i,j]), adj = if (j == 1) 0 else 0.5, cex = 0.9)
  }
  segments(0.04, 0.36, 0.97, 0.36, lwd = 1.2)
  text(0.04, 0.30, paste(strwrap(caption, width = 125), collapse = "\n"),
       adj = c(0, 1), cex = 0.8)
}
png(file.path(table_dir, "training_scenario_comparison_2026.png"),
    width = 10, height = 5.6, units = "in", res = 300)
draw_table(); dev.off()

# Forest plot: adjusted intervals are the primary comparison.
draw_figure <- function() {
  par(mar = c(6, 7, 3.5, 1.5), family = "serif", las = 1)
  y <- rev(seq_len(nrow(comparison)))
  limits <- range(c(0, comparison$DELTA_CI98_75_LOW, comparison$DELTA_CI98_75_HIGH))
  padding <- diff(limits) * 0.1
  plot(comparison$DELTA_R, y, type = "n", xlim = limits + c(-padding, padding),
       ylim = c(0.6, nrow(comparison) + 0.7), yaxt = "n", ylab = "", xlab = "",
       bty = "n", main = "Change in prediction ability with expanded training", cex.main = 1.05)
  axis(2, at = y, labels = unname(trait_labels[comparison$TRAIT]), tick = FALSE)
  abline(v = 0, lty = 2, col = "gray45")
  segments(comparison$DELTA_CI98_75_LOW, y, comparison$DELTA_CI98_75_HIGH, y,
           lwd = 2, col = "#24536B")
  segments(comparison$DELTA_CI98_75_LOW, y-0.07, comparison$DELTA_CI98_75_LOW, y+0.07, col="#24536B")
  segments(comparison$DELTA_CI98_75_HIGH, y-0.07, comparison$DELTA_CI98_75_HIGH, y+0.07, col="#24536B")
  points(comparison$DELTA_R, y, pch = 19, cex = 1.2, col = "#24536B")
  text(comparison$DELTA_R, y+0.23, sprintf("%+.3f", comparison$DELTA_R), cex = 0.85)
  mtext("Difference in Pearson correlation (expanded - historical)", side = 1, line = 2.4)
  mtext(paste0("Bonferroni-adjusted ", 100*(1-alpha/length(trait_order)),
               "% bootstrap intervals; ", format(n_boot,big.mark=","), " paired resamples"),
        side = 1, line = 4.1, cex = 0.8)
}
png(file.path(figure_dir, "training_scenario_comparison_2026.png"),
    width = 7.5, height = 4.8, units = "in", res = 300)
draw_figure(); dev.off()
print(thesis_table, row.names = FALSE)
message("Thesis table saved in: ", table_dir, "; figure saved in: ", figure_dir)
