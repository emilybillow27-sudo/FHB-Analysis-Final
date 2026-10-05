save_figure <- function(plot, name, width = 9, height = 7) {
  ggplot2::ggsave(result_path(file.path("figures", name)), plot,
                 width = width, height = height, dpi = 300, bg = "white")
}
theme_fhb <- function() ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

figure_analysis <- function() {
  markers <- load_result("markers.rds")
  historical <- load_result("historical_me.rds")
  testing <- load_result("testing_se.rds")
  evaluation <- load_result("evaluation_me.rds")
  forward <- load_result("forward_predictions.rds")
  # predict.mmes returns response-scale predictions, including the intercept.
  # Use the saved values directly; do not add the training mean again.
  plot <- ggplot2::ggplot(forward[forward$YEAR == 2025, ], ggplot2::aes(BLUE, GBLUP)) +
    ggplot2::geom_point(color = "#2166AC", alpha = 0.7) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::facet_wrap(~TRAIT, scales = "free") +
    ggplot2::labs(title = "2025 forward prediction", x = "Observed adjusted genotype mean", y = "Predicted phenotype") + theme_fhb()
  save_figure(plot, "forward_prediction_2025.png")
  plot <- ggplot2::ggplot(forward[forward$YEAR == 2026 & forward$environment == "Across locations", ],
                          ggplot2::aes(BLUE, GBLUP, color = TRAINING_POPULATION)) +
    ggplot2::geom_point(alpha = 0.7) + ggplot2::facet_wrap(~TRAIT, scales = "free") +
    ggplot2::labs(title = "2026 forward prediction", x = "Observed adjusted genotype mean", y = "Predicted phenotype", color = "Training population") + theme_fhb()
  save_figure(plot, "forward_prediction_2026.png")

  variable <- apply(markers, 2, stats::sd) > 0
  if (sum(variable) < 2L) stop("Too few variable markers for PCA.")
  # Only the first two scores are needed; retain all singular values for the scree plot.
  pca <- stats::prcomp(markers[, variable, drop = FALSE], center = TRUE, scale. = TRUE, rank. = 2L)
  variance <- pca$sdev^2 / sum(pca$sdev^2) * 100
  ids <- rownames(markers)
  in_train <- ids %in% historical$FullSampleName
  in_test <- ids %in% testing$FullSampleName
  set <- ifelse(in_train & in_test, "Training and testing", ifelse(in_train, "Training", ifelse(in_test, "Testing", "Other genotyped")))
  prefix <- substr(ids, 1, 2)
  prefix[!prefix %in% c("CO", "KS", "MT", "NE", "OK", "SD", "TX", "VA")] <- "Other"
  scores <- data.frame(FullSampleName = ids, PC1 = pca$x[, 1], PC2 = pca$x[, 2], program = prefix, population = set)
  write_table(scores, "PCA_scores.csv")
  scree <- data.frame(PC = seq_along(variance), Variance_explained = variance)
  write_table(scree, "PCA_variance_explained.csv")
  plot <- ggplot2::ggplot(scores, ggplot2::aes(PC1, PC2, color = program, shape = population)) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::labs(title = "Genetic structure of genotyped wheat lines", x = sprintf("PC1 (%.1f%%)", variance[1]),
                  y = sprintf("PC2 (%.1f%%)", variance[2]), color = "Program", shape = "Population") + theme_fhb()
  save_figure(plot, "PCA.png")
  plot <- ggplot2::ggplot(scree[scree$PC <= 64, ], ggplot2::aes(PC, Variance_explained)) +
    ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(x = "Principal component", y = "Variance explained (%)") + theme_fhb()
  save_figure(plot, "PCA_scree.png")

  correlations <- list()
  for (population in c("Historical", "2025 UIUC", "2026 across locations")) {
    dat <- switch(population, Historical = historical, `2025 UIUC` = testing, `2026 across locations` = evaluation)
    wide <- tidyr::pivot_wider(dplyr::summarise(dplyr::group_by(dat, genotype, TRAIT),
      value = mean(predicted.value), .groups = "drop"), names_from = TRAIT, values_from = value)
    traits <- intersect(config$traits, names(wide))
    m <- as.matrix(wide[traits])
    pairs <- expand.grid(Trait_1 = traits, Trait_2 = traits, stringsAsFactors = FALSE)
    pairs$Pearson_r <- mapply(function(a, b) safe_cor(m[, a], m[, b]), pairs$Trait_1, pairs$Trait_2)
    pairs$Paired_entries <- mapply(function(a, b) sum(is.finite(m[, a]) & is.finite(m[, b])), pairs$Trait_1, pairs$Trait_2)
    pairs$Population <- population
    correlations[[population]] <- pairs
  }
  correlations <- dplyr::bind_rows(correlations)
  write_table(correlations, "trait_correlations.csv")
  plot <- ggplot2::ggplot(correlations, ggplot2::aes(Trait_1, Trait_2, fill = Pearson_r)) +
    ggplot2::geom_tile(color = "white") + ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Pearson_r))) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1)) +
    ggplot2::facet_wrap(~Population) + ggplot2::coord_fixed() + ggplot2::labs(x = NULL, y = NULL, fill = "Pearson r") + theme_fhb()
  save_figure(plot, "trait_correlations.png", 13, 5)

  raw <- load_result("raw.rds")$testing
  observed <- dplyr::summarise(dplyr::group_by(long_phenotypes(raw), genotype, SUB_NURNAME, TRAIT),
                               value = mean(y, na.rm = TRUE), .groups = "drop")
  plot <- ggplot2::ggplot(observed[is.finite(observed$value), ], ggplot2::aes(value)) +
    ggplot2::geom_histogram(bins = 20, fill = "#2166AC", color = "white") +
    ggplot2::facet_wrap(~TRAIT, scales = "free") +
    ggplot2::labs(title = "2025 UIUC trait distributions", x = "Entry-nursery mean (DON: ppm; other traits: %)", y = "Entry-nursery count") + theme_fhb()
  save_figure(plot, "trait_distributions_2025.png")
  diagnostic_figures(markers)
}

# Diagnostics consolidated from the former GRM, LD, and variance-component scripts.
diagnostic_figures <- function(markers) {
  grm <- load_result("GRM.rds")
  cluster_order <- stats::hclust(stats::as.dist(1 - grm))$order
  ordered <- grm[cluster_order, cluster_order]
  cells <- expand.grid(row = seq_len(nrow(ordered)), column = seq_len(ncol(ordered)))
  cells$relationship <- as.vector(ordered)
  plot <- ggplot2::ggplot(cells, ggplot2::aes(row, column, fill = relationship)) +
    ggplot2::geom_raster() + ggplot2::coord_fixed() +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
    ggplot2::labs(title = "Clustered genomic relationships", x = "Genotypes", y = "Genotypes") +
    theme_fhb() + ggplot2::theme(axis.text = ggplot2::element_blank())
  save_figure(plot, "GRM_heatmap.png")

  vcf <- gaston::read.vcf(config$markers, convert.chr = FALSE)
  map <- vcf@snps[match(colnames(markers), vcf@snps$id), c("id", "chr", "pos")]
  if (anyNA(map)) stop("Marker positions do not match the PCA/GRM matrix.")
  map <- map[order(map$chr, map$pos), ]
  same <- as.character(map$chr[-1]) == as.character(map$chr[-nrow(map)])
  current <- map$id[-1][same]
  previous <- map$id[-nrow(map)][same]
  distance <- (map$pos[-1] - map$pos[-nrow(map)])[same]
  ld <- data.frame(dist_bp = distance, dist_kb = distance / 1000,
    r2 = mapply(function(a, b) safe_cor(markers[, a], markers[, b])^2, current, previous))
  write_table(ld, "ld_physical_adjacent.csv")
  plot <- ggplot2::ggplot(ld, ggplot2::aes(dist_kb, r2)) +
    ggplot2::geom_point(alpha = 0.25, size = 0.7) +
    ggplot2::labs(title = "Linkage disequilibrium between adjacent SNPs",
                  x = "Physical distance (kb)", y = expression(r^2)) + theme_fhb()
  save_figure(plot, "ld_physical_adjacent.png")

  h2 <- utils::read.csv(result_path("tables/heritabilities.csv"))
  variance <- tidyr::pivot_longer(h2, cols = c("Vg", "Vge", "Ve"), names_to = "Component", values_to = "Variance")
  plot <- ggplot2::ggplot(variance[is.finite(variance$Variance), ],
    ggplot2::aes(Location, Variance, fill = Component)) +
    ggplot2::geom_col(position = "dodge") + ggplot2::facet_grid(Population ~ Trait, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Estimated variance") + theme_fhb() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_figure(plot, "variance_components.png", 13, 9)
}
