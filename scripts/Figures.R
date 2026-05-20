# Load packages
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reshape2)
library(viridis)

# Define thesis plotting theme
theme_thesis <- theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey90"),
    plot.title = element_text(face = "bold")
  )

# Define traits
traits <- c("INC", "SEV", "DON")

# Create environment label for testing data
test <- test |>
  dplyr::mutate(Environment = paste(SUB_NURNAME, YR, sep = "_"))

# Plot raw phenotype histograms
raw_hist <- test |>
  tidyr::pivot_longer(cols = all_of(traits), names_to = "Trait", values_to = "Value") |>
  ggplot(aes(Value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  facet_wrap(~Trait, scales = "free") +
  labs(title = "Raw Phenotype Distributions", x = "Value", y = "Count") +
  theme_thesis

# Plot raw phenotype boxplots by environment
raw_box <- test |>
  tidyr::pivot_longer(cols = all_of(traits), names_to = "Trait", values_to = "Value") |>
  ggplot(aes(Environment, Value)) +
  geom_boxplot(outlier.alpha = 0.3) +
  facet_wrap(~Trait, scales = "free") +
  coord_flip() +
  labs(title = "Raw Phenotypes by Environment", x = "Environment", y = "Value") +
  theme_thesis

# Count observations per environment
env_counts <- test |>
  dplyr::group_by(Environment) |>
  dplyr::summarise(n = n(), .groups = "drop")

# Plot environment sample sizes
env_bar <- env_counts |>
  ggplot(aes(x = reorder(Environment, n), y = n)) +
  geom_col(fill = "darkolivegreen3") +
  coord_flip() +
  labs(title = "Sample Size per Environment", x = "Environment", y = "Number of Observations") +
  theme_thesis

# Compute raw phenotype correlation matrix
raw_corr <- test |>
  dplyr::select(all_of(traits)) |>
  cor(use = "pairwise.complete.obs")

# Plot raw phenotype correlation heatmap
corr_heatmap <- melt(raw_corr) |>
  ggplot(aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_viridis(option = "C", limits = c(-1, 1)) +
  labs(title = "Raw Phenotype Correlation Heatmap", x = "", y = "") +
  theme_thesis

# Plot within-environment BLUE histograms
blue_hist <- blues_se_train |>
  ggplot(aes(predicted.value)) +
  geom_histogram(bins = 30, fill = "tan3", color = "white") +
  facet_wrap(~TRAIT, scales = "free") +
  labs(title = "Within-Environment BLUE Distributions", x = "BLUE", y = "Count") +
  theme_thesis

# Plot standard error distributions for within-environment BLUEs
se_hist <- blues_se_train |>
  ggplot(aes(std.error)) +
  geom_histogram(bins = 30, fill = "orchid3", color = "white") +
  facet_wrap(~TRAIT, scales = "free") +
  labs(title = "Standard Error Distributions (Within-Environment BLUEs)", x = "Std. Error", y = "Count") +
  theme_thesis

# Prepare across-environment BLUEs for plotting
blue_pairs <- blues_me_train |>
  dplyr::group_by(FullSampleName, TRAIT) |>
  dplyr::summarise(predicted.value = mean(predicted.value, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = TRAIT, values_from = predicted.value)

# Plot across-environment BLUE histograms
me_hist <- blues_me_train |>
  ggplot(aes(predicted.value)) +
  geom_histogram(bins = 30, fill = "skyblue3", color = "white") +
  facet_wrap(~TRAIT, scales = "free") +
  labs(title = "Across-Environment BLUE Distributions", x = "BLUE", y = "Count") +
  theme_thesis

# Compute correlation matrix for across-environment BLUEs
blue_corr <- blue_pairs[, traits] |>
  cor(use = "pairwise.complete.obs")

# Plot across-environment BLUE correlation heatmap
blue_corr_heatmap <- melt(blue_corr) |>
  ggplot(aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_viridis(option = "C", limits = c(-1, 1)) +
  labs(title = "Across-Environment BLUE Correlation Heatmap", x = "", y = "") +
  theme_thesis

# Save figures
ggsave("raw_hist.png", raw_hist, width = 7, height = 5, dpi = 300)
ggsave("raw_box.png", raw_box, width = 7, height = 6, dpi = 300)
ggsave("env_bar.png", env_bar, width = 6, height = 6, dpi = 300)
ggsave("raw_corr_heatmap.png", corr_heatmap, width = 5, height = 4, dpi = 300)
ggsave("blue_hist.png", blue_hist, width = 7, height = 5, dpi = 300)
ggsave("se_hist.png", se_hist, width = 7, height = 5, dpi = 300)
ggsave("me_hist.png", me_hist, width = 7, height = 5, dpi = 300)
ggsave("blue_corr_heatmap.png", blue_corr_heatmap, width = 5, height = 4, dpi = 300)