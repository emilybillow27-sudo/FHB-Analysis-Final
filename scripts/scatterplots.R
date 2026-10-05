# Observed versus predicted scatterplots for 2025
# Author: Emily Billow
# Purpose: Compare observed 2025 UIUC adjusted genotype means with
# genomic estimated breeding values for INC, SEV, DON, and FDK

# Run from the project root (FHB Analysis Final).
if (!file.exists("data/FHB_Project_Training_Data.csv")) {
  stop("Set the working directory to the FHB Analysis Final project root.")
}
for (output_dir in c("results/intermediate", "results/tables", "results/figures")) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Read saved outputs from the main pipeline; no preloaded R objects are required.
forward_2025 <- read.csv("results/tables/forward_predictions_2025.csv", stringsAsFactors = FALSE)
blues_se_test <- readRDS("results/intermediate/blues_se_test.rds")
blues_me_train <- readRDS("results/intermediate/blues_me_train.rds")
traits <- c("INC", "SEV", "DON", "FDK")

# Confirm required columns
required_prediction_columns <- c(
  "FullSampleName",
  "GEBV",
  "TRAIT"
)

required_observed_columns <- c(
  "FullSampleName",
  "TRAIT",
  "predicted.value"
)

if (!all(required_prediction_columns %in% names(forward_2025))) {
  stop("forward_2025 is missing one or more required columns.")
}

if (!all(required_observed_columns %in% names(blues_se_test))) {
  stop("blues_se_test is missing one or more required columns.")
}

# Calculate one observed adjusted genotype mean for each genotype and trait
observed_2025 <- blues_se_test %>%
  dplyr::filter(TRAIT %in% traits) %>%
  dplyr::group_by(FullSampleName, TRAIT) %>%
  dplyr::summarise(
    observed = mean(predicted.value, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate historical training-population means
training_means <- blues_me_train %>%
  dplyr::filter(TRAIT %in% traits) %>%
  dplyr::group_by(TRAIT) %>%
  dplyr::summarise(
    trait_mean = mean(predicted.value, na.rm = TRUE),
    .groups = "drop"
  )

# Join predictions and observations
observed_predicted_2025 <- forward_2025 %>%
  dplyr::filter(TRAIT %in% traits) %>%
  dplyr::rename(
    predicted = GEBV
  ) %>%
  dplyr::inner_join(
    observed_2025,
    by = c("FullSampleName", "TRAIT")
  ) %>%
  dplyr::left_join(
    training_means,
    by = "TRAIT"
  ) %>%
  dplyr::filter(
    !is.na(observed),
    !is.na(predicted),
    !is.na(trait_mean)
  )

# Set this to TRUE only if GEBV is a centered genomic breeding value.
# If Zach confirms that GEBV already includes the population mean, use FALSE.
add_training_mean <- TRUE

if (add_training_mean) {
  observed_predicted_2025 <- observed_predicted_2025 %>%
    dplyr::mutate(
      predicted_on_observed_scale = predicted + trait_mean
    )
} else {
  observed_predicted_2025 <- observed_predicted_2025 %>%
    dplyr::mutate(
      predicted_on_observed_scale = predicted
    )
}

# Check the number of genotypes used for each trait
prediction_counts <- observed_predicted_2025 %>%
  dplyr::count(TRAIT, name = "n_genotypes")

print(prediction_counts)

# Save the joined prediction dataset
dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  observed_predicted_2025,
  "results/tables/observed_predicted_2025.csv",
  row.names = FALSE
)

# Create one scatterplot for a single trait
plot_observed_predicted <- function(dat, trait_name) {
  
  trait_data <- dat %>%
    dplyr::filter(TRAIT == trait_name)
  
  predictive_correlation <- cor(
    trait_data$observed,
    trait_data$predicted_on_observed_scale,
    use = "complete.obs"
  )
  
  rmse <- sqrt(
    mean(
      (
        trait_data$observed -
          trait_data$predicted_on_observed_scale
      )^2,
      na.rm = TRUE
    )
  )
  
  ggplot2::ggplot(
    trait_data,
    ggplot2::aes(
      x = observed,
      y = predicted_on_observed_scale
    )
  ) +
    ggplot2::geom_point(
      color = "#2166AC",
      size = 2.5,
      alpha = 0.75
    ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "black"
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      se = FALSE,
      color = "#B2182B"
    ) +
    ggplot2::labs(
      title = paste("2025 UIUC:", trait_name),
      subtitle = paste0(
        "Predictive correlation = ",
        round(predictive_correlation, 3),
        " | RMSE = ",
        round(rmse, 3),
        " | n = ",
        nrow(trait_data)
      ),
      x = "Observed adjusted genotype mean",
      y = "Predicted value on observed scale"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# Create output folder
dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# Generate and save individual trait plots
plots_2025 <- list()

for (trait_name in traits) {
  
  plots_2025[[trait_name]] <- plot_observed_predicted(
    observed_predicted_2025,
    trait_name
  )
  
  ggplot2::ggsave(
    filename = paste0(
      "results/figures/observed_vs_predicted_2025_",
      trait_name,
      ".png"
    ),
    plot = plots_2025[[trait_name]],
    width = 7,
    height = 6,
    dpi = 300
  )
}

# Combine all traits into one figure
combined_2025_scatterplots <- patchwork::wrap_plots(
  plots_2025,
  ncol = 2
)

# Display the combined figure
print(combined_2025_scatterplots)

# Save the combined figure
ggplot2::ggsave(
  "results/figures/observed_vs_predicted_2025_all_traits.png",
  combined_2025_scatterplots,
  width = 14,
  height = 11,
  dpi = 300
)