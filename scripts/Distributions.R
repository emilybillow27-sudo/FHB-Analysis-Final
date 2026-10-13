# Observed trait distributions: historical, 2025, and 2026
# Run from the project root:
# source("scripts/Distributions.R")
# Matches the taller 2025 UIUC figure; no prediction models are rerun.

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

project_dir <- getOption("fhb.project_dir", "/Users/emilybillow/Desktop/FHB Analysis Final")
figure_dir <- file.path(project_dir, "results", "figures", "distributions")
table_dir <- file.path(project_dir, "results", "tables", "distributions")
for (directory in c(figure_dir, table_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

trait_order <- c("INC", "SEV", "DON", "FDK")
trait_labels <- c(INC = "Incidence", SEV = "Severity", DON = "DON", FDK = "FDK")
trait_units <- c(INC = "%", SEV = "%", DON = "ppm", FDK = "%")
trait_colors <- c(INC = "#214A68", SEV = "#4E7FA3", DON = "#687985", FDK = "#7699B2")

average_observations <- function(raw) {
  if (anyNA(raw$ID) || anyNA(raw$Environment)) {
    stop("Missing genotype or environment identifiers.")
  }

  for (trait in setdiff(trait_order, names(raw))) {
    raw[[trait]] <- NA_real_
  }

  if (!all(vapply(raw[trait_order], is.numeric, logical(1)))) {
    stop("Trait columns must be numeric.")
  }

  raw %>%
    pivot_longer(
      all_of(trait_order),
      names_to = "TRAIT",
      values_to = "value"
    ) %>%
    group_by(Environment, ID, TRAIT) %>%
    summarise(
      n_observations = sum(is.finite(value)),
      value = if (any(is.finite(value))) {
        mean(value[is.finite(value)])
      } else {
        NA_real_
      },
      .groups = "drop"
    )
}

plot_distribution <- function(dat, trait, adjusted = FALSE) {
  vals <- dat %>% filter(TRAIT == trait, is.finite(value))
  n_missing <- sum(dat$TRAIT == trait & !is.finite(dat$value))
  if (!nrow(vals)) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Data not yet available",
                 size = 6, color = "#687985") +
        xlim(0, 1) + ylim(0, 1) +
        labs(title = trait_labels[[trait]], subtitle = "No observed measurements") +
        theme_void(base_size = 16) +
        theme(plot.title = element_text(face = "bold", color = "#214A68", size = 19),
              plot.subtitle = element_text(color = "#687985", size = 12),
              plot.margin = margin(12, 16, 12, 12))
    )
  }
  x_label <- paste0(if (adjusted) "BLUE" else "Mean observed value", " (", trait_units[[trait]], ")")
  # Fixed, interpretable bins per trait; DON retains its own ppm scale.
  bin_width <- if (trait == "DON") 5 else 10
  p <- ggplot(vals, aes(x = value)) +
    geom_histogram(binwidth = bin_width, boundary = 0, closed = "left",
                   fill = trait_colors[[trait]], color = "white", linewidth = 0.5) +
    labs(title = trait_labels[[trait]],
         subtitle = paste0("n = ", nrow(vals), " entry-environment means",
                           if (n_missing) paste0("; ", n_missing, " missing") else ""),
         x = x_label, y = "Genotype-environment observations") +
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


save_distribution <- function(dat, stem, title, subtitle, caption) {
  stopifnot(!anyDuplicated(dat[c("ID", "Environment", "TRAIT")]))

  plots <- lapply(trait_order, function(trait) {
    plot_distribution(dat, trait)
  })

  fig <- wrap_plots(plots, ncol = 2) +
    plot_annotation(
      title = title,
      subtitle = subtitle,
      caption = caption,
      theme = theme(
        plot.title = element_text(face = "bold", size = 24, color = "#214A68"),
        plot.subtitle = element_text(size = 16, color = "#687985"),
        plot.caption = element_text(size = 11, hjust = 0, color = "#687985"),
        plot.background = element_rect(fill = "white", color = NA)
      )
    )

  ggsave(
    file.path(figure_dir, paste0(stem, ".png")),
    fig,
    width = 10.5,
    height = 9,
    dpi = 300,
    bg = "white"
  )

  write.csv(
    dat,
    file.path(table_dir, paste0(stem, "_values.csv")),
    row.names = FALSE,
    na = ""
  )

  message("Saved ", stem)
  invisible(fig)
}

# Match the 2025 descriptive figure: all five UIUC sub-nurseries.
# Entries shared across sub-nurseries remain separate.

uiuc_2025_raw <- read.csv(
  file.path(project_dir, "data", "FHB_Project_Testing_Data_2025.csv"),
  stringsAsFactors = FALSE
)

uiuc_2025_raw <- uiuc_2025_raw %>%
  filter(toupper(trimws(EXPT)) == "UIUC", YR == 2025) %>%
  mutate(Environment = SUB_NURNAME)

uiuc_2025_means <- average_observations(uiuc_2025_raw)

save_distribution(
  uiuc_2025_means,
  "UIUC_2025_observed_distributions_taller",
  "2025 UIUC trait distributions",
  "Observed means within entry and sub-nursery",
  paste0(
    "Replicates averaged within entry and sub-nursery at UIUC in 2025 (",
    n_distinct(uiuc_2025_raw$Environment), " sub-nurseries).\n",
    "Entries shared across sub-nurseries remain separate. Units and axis ranges vary."
  )
)

historical_raw <- read.csv(
  file.path(project_dir, "data", "FHB_Project_Training_Data.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

historical <- historical_raw %>%
  transmute(
    ID = germplasmName,
    Environment = studyName,
    INC = `FHB.incidence.....CO_321.0001149`,
    SEV = `FHB.severity.....CO_321.0001440`,
    DON = `FHB.DON.content...ppm.CO_321.0001154`,
    FDK = `FHB.grain.incidence.....CO_321.0001155`
  )

historical_means <- average_observations(historical)

save_distribution(
  historical_means,
  "Historical_observed_distributions",
  "Historical training trait distributions",
  "Observed means within genotype and study",
  paste0(
    "Replicates averaged within genotype and study; ",
    n_distinct(historical$Environment), " historical studies.\n",
    "Genotypes shared across studies remain separate. Units and axis ranges vary."
  )
)

# Use all wheat entries, matching the descriptive scope of the 2025 figure.
# Exclude Triticale; do not restrict to the 60 forward-validation genotypes.
for (location in c("UNL", "UIUC", "SDSU")) {
  raw <- read.csv(
    file.path(project_dir, "data", paste0("FHB_Project_Testing_Data_", location, "_2026.csv")),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  raw <- raw[!is.na(raw$ID) & raw$ID != "Triticale", ]
  raw$Environment <- location
  raw$INC <- raw$Incidence
  raw$SEV <- raw$Severity
  observed <- average_observations(raw)

  save_distribution(
    observed,
    paste0(location, "_2026_observed_distributions"),
    paste0("2026 ", location, " trait distributions"),
    "Observed means within wheat entry and location",
    paste0(
      "Replicates averaged within entry at ", location, "; all available wheat entries.\n",
      "Triticale excluded. Unavailable traits remain missing. Units and axis ranges vary."
    )
  )
}
