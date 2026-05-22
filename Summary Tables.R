library(tidyverse)
library(e1071)

# ============================================================
# FUNCTION: Phenotype Summary Table (mean, SD, range, skew)
# ============================================================

phenotype_summary <- function(df, traits) {
  
  df_long <- df %>%
    pivot_longer(
      cols = all_of(traits),
      names_to = "trait",
      values_to = "value"
    )
  
  df_long %>%
    group_by(trait) %>%
    summarise(
      n = sum(!is.na(value)),
      mean = mean(value, na.rm = TRUE),
      sd   = sd(value, na.rm = TRUE),
      min  = min(value, na.rm = TRUE),
      max  = max(value, na.rm = TRUE),
      range = max - min,
      skewness = skewness(value, na.rm = TRUE, type = 2),
      .groups = "drop"
    )
}

# ============================================================
# TESTING POPULATION SUMMARY
# ============================================================

testing <- read_csv("data/FHB_Project_Testing_Data.csv")

testing_summary <- phenotype_summary(
  testing,
  traits = c("INC", "SEV", "DON")
) %>%
  mutate(dataset = "Testing")

# ============================================================
# TRAINING POPULATION SUMMARY
# ============================================================

train <- read_csv("data/FHB_Project_Training_Data.csv") %>%
  rename(
    INC = `FHB.incidence.....CO_321.0001149`,
    SEV = `FHB.severity.....CO_321.0001440`,
    DON = `FHB.DON.content...ppm.CO_321.0001154`,
    DI  = `FHB.disease.index.....CO_321.0501030`,
    ID  = uID
  )

training_summary <- phenotype_summary(
  train,
  traits = c("INC", "SEV", "DON", "DI")
) %>%
  mutate(dataset = "Training")

# ============================================================
# MERGE WITH VISUAL SEPARATOR
# ============================================================

phenotype_summary_full <- bind_rows(
  testing_summary,
  separator_row,
  training_summary
)

# ============================================================
# SAVE TABLE
# ============================================================

write_csv(
  phenotype_summary_full,
  "results/tables/phenotype_summary_table.csv"
)

phenotype_summary_full