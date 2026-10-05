# Combined 2025 UIUC trait distributions
# Run independently; the full genomic prediction pipeline need not be rerun.
# Packages: dplyr, tidyr, ggplot2, patchwork (already available locally).
# Example from the project root:
# source("/path/to/plot_uiuc_2025_distributions.R")
# Set these options before sourcing to use other locations:
# options(fhb.project_dir = "/path/to/FHB Analysis Final",
#         fhb.figure_dir = "/path/to/figures")

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

project_dir <- getOption("fhb.project_dir", "/Users/emilybillow/Desktop/FHB Analysis Final")
figure_dir <- getOption("fhb.figure_dir", file.path(project_dir, "results", "figures"))
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
trait_order <- c("INC", "SEV", "DON", "FDK")
trait_labels <- c(INC = "Incidence", SEV = "Severity", DON = "Deoxynivalenol (DON)",
                  FDK = "Fusarium-damaged kernels")
trait_units <- c(INC = "%", SEV = "%", DON = "ppm", FDK = "%")
trait_colors <- c(INC = "#214A68", SEV = "#4E7FA3", DON = "#687985", FDK = "#7699B2")

raw <- read.csv(file.path(project_dir, "data", "FHB_Project_Testing_Data_2025.csv"),
                stringsAsFactors = FALSE)
stopifnot(all(c("ID", "SUB_NURNAME", "EXPT", "YR", trait_order) %in% names(raw)))
raw <- raw %>% filter(toupper(trimws(EXPT)) == "UIUC", YR == 2025)
if (!nrow(raw)) stop("No 2025 UIUC observations found.")
if (any(is.na(raw$ID) | trimws(as.character(raw$ID)) == "")) stop("Missing entry ID.")
if (any(is.na(raw$SUB_NURNAME) | trimws(raw$SUB_NURNAME) == "")) stop("Missing nursery ID.")
if (!all(vapply(raw[trait_order], is.numeric, logical(1)))) stop("Trait columns must be numeric.")

# Use ID, not FullSampleName: Topcross entries have missing FullSampleName.
# Average replicates within ID and nursery. Keep the same entry in different
# nurseries separate; the unit is an entry-nursery, not a unique genotype.
# Unreplicated entries remain as their one observed value. Do not impute missing traits.
observed <- raw %>%
  pivot_longer(all_of(trait_order), names_to = "TRAIT", values_to = "value") %>%
  group_by(SUB_NURNAME, ID, TRAIT) %>%
  summarise(n_observations = sum(is.finite(value)),
            value = if (any(is.finite(value))) mean(value[is.finite(value)]) else NA_real_,
            .groups = "drop")

plot_distribution <- function(dat, trait, adjusted = FALSE) {
  vals <- dat %>% filter(TRAIT == trait, is.finite(value))
  n_missing <- sum(dat$TRAIT == trait & !is.finite(dat$value))
  x_label <- paste0(if (adjusted) "BLUE" else "Mean observed value", " (", trait_units[[trait]], ")")
  # Fixed, interpretable bins per trait; DON retains its own ppm scale.
  bin_width <- if (trait == "DON") 5 else 10
  p <- ggplot(vals, aes(x = value)) +
    geom_histogram(binwidth = bin_width, boundary = 0, closed = "left",
                   fill = trait_colors[[trait]], color = "white", linewidth = 0.5) +
    labs(title = trait_labels[[trait]],
         subtitle = paste0("n = ", nrow(vals), " entry-nursery values",
                           if (n_missing) paste0("; ", n_missing, " missing") else ""),
         x = x_label, y = "Number of entries") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    theme_minimal(base_size = 16, base_family = "sans") +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          plot.title = element_text(face = "bold", color = "#214A68", size = 19),
          plot.subtitle = element_text(color = "#687985", size = 12),
          axis.title = element_text(size = 14), axis.text = element_text(color = "#334652"),
          plot.margin = margin(12, 16, 12, 12))
  if (!nrow(vals)) p <- p + annotate("text", x = 0, y = 0, label = "No available values")
  p
}

save_distributions <- function(dat, stem, adjusted = FALSE) {
  if (!all(c("ID", "SUB_NURNAME", "TRAIT", "value") %in% names(dat))) stop("Missing plot columns.")
  if (anyDuplicated(dat[c("ID", "SUB_NURNAME", "TRAIT")])) stop("Duplicate entry-nursery-trait estimates.")
  plots <- lapply(trait_order, function(trait) plot_distribution(dat, trait, adjusted))
  caption <- if (adjusted) {
    "Within-nursery BLUEs from replicated Elite and WSS AYN trials. Unreplicated nurseries excluded.\nNegative adjusted estimates are retained. Axes use trait-specific units and ranges."
  } else {
    "Replicates averaged within each entry and nursery; all five nurseries included.\nRepeated entries in different nurseries remain separate. Axes use trait-specific units and ranges."
  }
  fig <- wrap_plots(plots, ncol = 2) +
    plot_annotation(title = "Trait distributions in the 2025 UIUC population",
                    subtitle = if (adjusted) "Adjusted entry estimates from replicated nurseries" else "Observed entry means across nurseries",
                    caption = caption,
                    theme = theme(plot.title = element_text(face = "bold", size = 24, color = "#214A68"),
                                  plot.subtitle = element_text(size = 16, color = "#687985"),
                                  plot.caption = element_text(size = 11, hjust = 0, color = "#687985"),
                                  plot.background = element_rect(fill = "white", color = NA)))
  ggsave(file.path(figure_dir, paste0(stem, ".png")), fig,
         width = 13, height = 8.5, units = "in", dpi = 300, bg = "white")
  write.csv(dat, file.path(figure_dir, paste0(stem, "_values.csv")), row.names = FALSE)
  print(dat %>% group_by(TRAIT) %>% summarise(n = sum(is.finite(value)),
        missing = sum(!is.finite(value)), .groups = "drop"))
  invisible(fig)
}

uiuc_2025_observed_figure <- save_distributions(observed, "UIUC_2025_observed_distributions")

# The BLUE version uses your existing saved estimates; no model is refitted.
blue_path <- file.path(project_dir, "results", "blues_se_test.rds")
if (file.exists(blue_path)) {
  blue <- readRDS(blue_path)
  stopifnot(all(c("ID", "SUB_NURNAME", "TRAIT", "predicted.value", "adjusted") %in% names(blue)))
  if (!all(blue$SUB_NURNAME %in% c("Elite", "WSS AYN"))) stop("Check BLUE nursery membership before plotting.")
  adjusted_values <- blue %>% filter(TRAIT %in% trait_order, adjusted == TRUE) %>%
    transmute(SUB_NURNAME, ID, TRAIT, value = predicted.value)
  uiuc_2025_blue_figure <- save_distributions(adjusted_values, "UIUC_2025_BLUE_distributions", TRUE)
} else {
  message("BLUE file absent; observed distributions saved. Run the 2025 BLUE section first for the adjusted version.")
}
message("Figures saved in: ", normalizePath(figure_dir))
