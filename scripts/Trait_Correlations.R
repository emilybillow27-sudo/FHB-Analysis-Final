# ============================================================
# Trait Correlation Heatmaps for INC, SEV, DON
# Using:
#   - blues_se_test  (testing BLUEs, single-environment)
#   - blues_me_train (training BLUEs, across-environment)
# ============================================================

library(tidyverse)
library(reshape2)
library(patchwork)

traits <- c("INC", "SEV", "DON")

# ------------------------------------------------------------
# 1. CLEAN + SUMMARIZE TESTING BLUEs
# ------------------------------------------------------------
blues_test_wide <- blues_se_test |>
  group_by(FullSampleName, TRAIT) |>
  summarise(pred = mean(predicted.value, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(
    names_from = TRAIT,
    values_from = pred
  )

# ------------------------------------------------------------
# 2. CLEAN + SUMMARIZE TRAINING BLUEs (ACROSS-ENV BLUEs)
# ------------------------------------------------------------
blues_train_wide <- blues_me_train |>
  group_by(FullSampleName, TRAIT) |>
  summarise(pred = mean(predicted.value, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(
    names_from = TRAIT,
    values_from = pred
  )

# ------------------------------------------------------------
# 3. Compute correlation matrices
# ------------------------------------------------------------
test_numeric  <- blues_test_wide  |> dplyr::select(all_of(traits))
train_numeric <- blues_train_wide |> dplyr::select(all_of(traits))

cor_test  <- cor(test_numeric,  use = "pairwise.complete.obs")
cor_train <- cor(train_numeric, use = "pairwise.complete.obs")

# ------------------------------------------------------------
# 4. Heatmap plotting function (red/blue + labels)
# ------------------------------------------------------------
plot_corr_heatmap <- function(cor_matrix, title_text) {
  
  melt(cor_matrix) |>
    ggplot(aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", value)), size = 5) +
    scale_fill_gradient2(
      low = "#2166ac",   # blue
      mid = "white",
      high = "#b2182b",  # red
      midpoint = 0,
      limits = c(-1, 1),
      name = "Correlation"
    ) +
    labs(title = title_text, x = "", y = "") +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
}

# ------------------------------------------------------------
# 5. Generate both heatmaps
# ------------------------------------------------------------
p_test <- plot_corr_heatmap(
  cor_test,
  "Correlation Heatmap for FHB Traits (Testing BLUEs)"
)

p_train <- plot_corr_heatmap(
  cor_train,
  "Correlation Heatmap for FHB Traits (Training BLUEs)"
)

# ------------------------------------------------------------
# 6. Combine vertically for thesis figure
# ------------------------------------------------------------
combined_corr <- p_test / p_train

# ------------------------------------------------------------
# 7. Display in RStudio
# ------------------------------------------------------------
combined_corr
