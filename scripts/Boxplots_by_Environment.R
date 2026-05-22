# ============================================================
# Boxplots by Environment for Training + Testing Populations
# ============================================================

library(tidyverse)

# ------------------------------------------------------------
# 0. Load phenotype data
# ------------------------------------------------------------

fhb   <- read_csv("data/FHB_Project_Testing_Data.csv")
train <- read_csv("data/FHB_Project_Training_Data.csv")

# ------------------------------------------------------------
# 1. Clean environment variables
# ------------------------------------------------------------

fhb <- fhb %>%
  mutate(EXPT = factor(EXPT, levels = unique(EXPT)))

train <- train %>%
  mutate(locationName = factor(locationName, levels = unique(locationName)))

# ------------------------------------------------------------
# 2. RENAME TESTING TRAIT COLUMNS
# ------------------------------------------------------------

# Testing dataset already has INC, SEV, DON — ensure they are correctly named
don_col_fhb <- grep("DON", colnames(fhb), value = TRUE)[1]
inc_col_fhb <- grep("INC", colnames(fhb), value = TRUE)[1]
sev_col_fhb <- grep("SEV", colnames(fhb), value = TRUE)[1]

names(fhb)[names(fhb) == don_col_fhb] <- "DON"
names(fhb)[names(fhb) == inc_col_fhb] <- "INC"
names(fhb)[names(fhb) == sev_col_fhb] <- "SEV"

fhb$ID <- as.character(fhb$ID)

# ------------------------------------------------------------
# 3. RENAME TRAINING TRAIT COLUMNS (ontology → INC/SEV/DON)
# ------------------------------------------------------------

don_col <- grep("DON.content", colnames(train), value = TRUE)[1]
inc_col <- grep("grain.incidence", colnames(train), value = TRUE)[1]
inc2_col <- grep("FHB.incidence", colnames(train), value = TRUE)[1]
sev_col <- grep("severity", colnames(train), value = TRUE)[1]
di_col  <- grep("disease.index", colnames(train), value = TRUE)[1]

# Use grain incidence if available, otherwise FHB.incidence
if (!is.na(inc_col)) {
  names(train)[names(train) == inc_col] <- "INC"
} else {
  names(train)[names(train) == inc2_col] <- "INC"
}

names(train)[names(train) == don_col] <- "DON"
names(train)[names(train) == sev_col] <- "SEV"
names(train)[names(train) == di_col]  <- "DI"

names(train)[names(train) == "uID"] <- "ID"
train$ID <- as.character(train$ID)

# ------------------------------------------------------------
# 4. Convert both datasets to long format
# ------------------------------------------------------------

fhb_long <- fhb %>%
  dplyr::select(EXPT, INC, SEV, DON) %>%
  pivot_longer(
    cols = c(INC, SEV, DON),
    names_to = "Trait",
    values_to = "Value"
  ) %>%
  mutate(Pop = "Testing")

train_long <- train %>%
  dplyr::select(locationName, INC, SEV, DON) %>%
  pivot_longer(
    cols = c(INC, SEV, DON),
    names_to = "Trait",
    values_to = "Value"
  ) %>%
  rename(EXPT = locationName) %>%
  mutate(Pop = "Training")

# ------------------------------------------------------------
# 5. Combine and remove NA rows (prevents ggplot warnings)
# ------------------------------------------------------------

combined_long <- bind_rows(fhb_long, train_long) %>%
  filter(!is.na(Value))

# ------------------------------------------------------------
# 6. Boxplots by environment
# ------------------------------------------------------------

p_box_env <- ggplot(combined_long, aes(x = EXPT, y = Value, fill = EXPT)) +
  geom_boxplot(outlier.color = "red", outlier.alpha = 0.7) +
  facet_grid(Pop ~ Trait, scales = "free_y") +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "Trait Distributions by Environment and Population",
    x = "Environment",
    y = "Trait Value"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ------------------------------------------------------------
# 7. Save figure
# ------------------------------------------------------------

ggsave("results/figures/boxplots_by_environment_training_testing.png",
       p_box_env, width = 12, height = 8, dpi = 300)
