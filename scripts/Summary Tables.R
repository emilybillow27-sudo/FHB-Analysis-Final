# Summary_Tables.R
# Author: Emily Billow
# Purpose: Summarize FHB phenotypes across the historical training population,
# the 2025 UIUC testing population, and the 2026 UNL and UIUC populations.

library(dplyr)
library(tidyr)
library(readr)
library(e1071)

# Traits expected in the testing datasets

testing_traits <- c(
  "INC",
  "SEV",
  "DON",
  "FDK"
)

# Traits expected in the historical dataset

historical_traits <- c(
  "INC",
  "SEV",
  "DON",
  "FDK",
  "DI"
)

# Calculate summary statistics for all available traits

phenotype_summary <- function(
    df,
    dataset_name,
    expected_traits
) {
  
  available_traits <- intersect(
    expected_traits,
    names(df)
  )
  
  if (length(available_traits) == 0) {
    stop(
      "No expected trait columns were found in ",
      dataset_name
    )
  }
  
  missing_traits <- setdiff(
    expected_traits,
    available_traits
  )
  
  if (length(missing_traits) > 0) {
    message(
      dataset_name,
      " does not contain: ",
      paste(missing_traits, collapse = ", ")
    )
  }
  
  df %>%
    dplyr::select(
      dplyr::all_of(available_traits)
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "Trait",
      values_to = "Value"
    ) %>%
    dplyr::group_by(Trait) %>%
    dplyr::summarise(
      n = sum(!is.na(Value)),
      mean = if (all(is.na(Value))) {
        NA_real_
      } else {
        mean(Value, na.rm = TRUE)
      },
      sd = if (sum(!is.na(Value)) < 2) {
        NA_real_
      } else {
        sd(Value, na.rm = TRUE)
      },
      minimum = if (all(is.na(Value))) {
        NA_real_
      } else {
        min(Value, na.rm = TRUE)
      },
      maximum = if (all(is.na(Value))) {
        NA_real_
      } else {
        max(Value, na.rm = TRUE)
      },
      skewness = if (sum(!is.na(Value)) < 3) {
        NA_real_
      } else {
        e1071::skewness(
          Value,
          na.rm = TRUE,
          type = 2
        )
      },
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      range = maximum - minimum,
      Dataset = dataset_name
    ) %>%
    dplyr::select(
      Dataset,
      Trait,
      n,
      mean,
      sd,
      minimum,
      maximum,
      range,
      skewness
    )
}

# Load the 2025 UIUC testing dataset

testing_2025 <- readr::read_csv(
  "data/FHB_Project_Testing_Data_2025.csv",
  show_col_types = FALSE
)

# Summarize the 2025 UIUC data

summary_2025 <- phenotype_summary(
  df = testing_2025,
  dataset_name = "2025 UIUC testing population",
  expected_traits = testing_traits
)

# Load the 2026 UNL testing dataset

testing_2026_unl <- readr::read_csv(
  "data/FHB_Project_Testing_Data_UNL_2026.csv",
  show_col_types = FALSE
) %>%
  dplyr::rename(
    INC = Incidence,
    SEV = Severity
  )

# Summarize the 2026 UNL data

summary_2026_unl <- phenotype_summary(
  df = testing_2026_unl,
  dataset_name = "2026 UNL testing population",
  expected_traits = testing_traits
)

# Load the 2026 UIUC testing dataset

testing_2026_uiuc <- readr::read_csv(
  "data/FHB_Project_Testing_Data_UIUC_2026.csv",
  show_col_types = FALSE
) %>%
  dplyr::rename(
    INC = Incidence,
    SEV = Severity
  )

# Summarize the 2026 UIUC data

summary_2026_uiuc <- phenotype_summary(
  df = testing_2026_uiuc,
  dataset_name = "2026 UIUC testing population",
  expected_traits = testing_traits
)

# Load the historical training dataset

training <- readr::read_csv(
  "data/FHB_Project_Training_Data.csv",
  show_col_types = FALSE
)

# Rename historical phenotype columns

historical_rename <- c(
  "FHB.incidence.....CO_321.0001149" =
    "INC",
  "FHB.severity.....CO_321.0001440" =
    "SEV",
  "FHB.DON.content...ppm.CO_321.0001154" =
    "DON",
  "FHB.grain.incidence.....CO_321.0001155" =
    "FDK",
  "FHB.disease.index.....CO_321.0501030" =
    "DI",
  "uID" =
    "ID"
)

for (old_name in names(historical_rename)) {
  
  new_name <- historical_rename[[old_name]]
  
  if (old_name %in% names(training)) {
    names(training)[
      names(training) == old_name
    ] <- new_name
  }
}

# Summarize the historical training data

summary_historical <- phenotype_summary(
  df = training,
  dataset_name = "Historical training population",
  expected_traits = historical_traits
)

# Combine all population summaries

phenotype_summary_full <- dplyr::bind_rows(
  summary_historical,
  summary_2025,
  summary_2026_unl,
  summary_2026_uiuc
)

# Set a consistent dataset order

phenotype_summary_full <- phenotype_summary_full %>%
  dplyr::mutate(
    Dataset = factor(
      Dataset,
      levels = c(
        "Historical training population",
        "2025 UIUC testing population",
        "2026 UNL testing population",
        "2026 UIUC testing population"
      )
    ),
    Trait = factor(
      Trait,
      levels = c("INC", "SEV", "DON", "FDK", "DI")
    )
  ) %>%
  dplyr::arrange(
    Dataset,
    Trait
  ) %>%
  dplyr::mutate(
    Dataset = as.character(Dataset),
    Trait = as.character(Trait)
  )

# Create output directories

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

# Save the combined summary table

readr::write_csv(
  phenotype_summary_full,
  "results/tables/phenotype_summary_all_populations.csv"
)

# Print the complete table

print(phenotype_summary_full)

# Print each dataset separately for easier inspection

cat("\nHistorical training population:\n")
print(summary_historical)

cat("\n2025 UIUC testing population:\n")
print(summary_2025)

cat("\n2026 UNL testing population:\n")
print(summary_2026_unl)

cat("\n2026 UIUC testing population:\n")
print(summary_2026_uiuc)