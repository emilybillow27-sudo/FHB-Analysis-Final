# ============================================================
# Testing Population Heritabilities
# ============================================================

library(tidyverse)
library(lme4)

# Import dataset
fhb <- read_csv("FHB_Project_Testing_Data.csv")

# Function for testing population (always 2 reps)
calc_h2_test <- function(trait) {
  df <- fhb %>% filter(!is.na(.data[[trait]]))
  n_rep <- n_distinct(df$REP)
  n_obs <- nrow(df)

  # Fit model with genotype + replicate
  model <- lmer(as.formula(paste0(trait, " ~ (1|ID) + (1|REP)")), data = df)

  varcomp <- as.data.frame(VarCorr(model))
  Vg <- varcomp$vcov[varcomp$grp == "ID"]
  Ve <- varcomp$vcov[varcomp$grp == "Residual"]

  # Broad-sense H2 for replicated trial
  H2 <- Vg / (Vg + Ve / n_rep)

  # Standard error of H2
  SE <- sqrt(2 * (H2^2) / n_obs)

  tibble(
    Trait = trait,
    H2 = H2,
    SE = SE,
    Reps = n_rep
  )
}

traits <- c("DON", "INC", "SEV")
testing_h2 <- map_df(traits, calc_h2_test)
testing_h2


# Plot Testing H2 with SE bars
ggplot(testing_h2, aes(x = Trait, y = H2, fill = Trait)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = H2 - SE, ymax = H2 + SE),
    width = 0.15,
    linewidth = 0.8
  ) +
  scale_fill_manual(values = c(
    DON = "#6A8F72",
    INC = "#D6B36A",
    SEV = "#C9C1A3"
  )) +
  labs(
    title = "Testing Population Heritability (H²)",
    y = "Heritability (H²)",
    x = NULL
  ) +
  theme_minimal(base_size = 14)


# ============================================================
# Training Population Heritabilities
# ============================================================

train <- read_csv("FHB_Project_Training_Data.csv")

train <- train %>%
  rename(
    DON = `FHB.DON.content...ppm.CO_321.0001154`,
    INC = `FHB.incidence.....CO_321.0001149`,
    SEV = `FHB.severity.....CO_321.0001440`,
    ID = uID
  )

# Function that adapts to number of reps per trait
calc_h2_train <- function(trait) {
  df <- train %>% filter(!is.na(.data[[trait]]))
  n_rep <- n_distinct(df$replicate)
  n_obs <- nrow(df)

  # Case 1: replicated trial (≥2 reps)
  if (n_rep >= 2) {
    model <- lmer(as.formula(paste0(trait, " ~ (1|ID) + (1|replicate)")), data = df)

    vc <- as.data.frame(VarCorr(model))
    Vg <- vc$vcov[vc$grp == "ID"]
    Ve <- vc$vcov[vc$grp == "Residual"]

    H2 <- Vg / (Vg + Ve / n_rep)
    SE <- sqrt(2 * (H2^2) / n_obs)

    return(tibble(Trait = trait, H2 = H2, SE = SE, Reps = n_rep))
  }

  # Case 2: unreplicated trial (1 rep)
  if (n_rep == 1) {
    model <- lmer(as.formula(paste0(trait, " ~ (1|ID)")), data = df)

    vc <- as.data.frame(VarCorr(model))
    Vg <- vc$vcov[vc$grp == "ID"]
    Ve <- vc$vcov[vc$grp == "Residual"]

    H2 <- Vg / (Vg + Ve) # repeatability
    SE <- sqrt(2 * (H2^2) / n_obs)

    return(tibble(Trait = trait, H2 = H2, SE = SE, Reps = n_rep))
  }
}

traits <- c("DON", "INC", "SEV")
training_h2 <- map_df(traits, calc_h2_train)
training_h2


# Plot Training H2 with SE bars
ggplot(training_h2, aes(x = Trait, y = H2, fill = Trait)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = H2 - SE, ymax = H2 + SE),
    width = 0.15,
    linewidth = 0.8
  ) +
  scale_fill_manual(values = c(
    DON = "#6A8F72",
    INC = "#D6B36A",
    SEV = "#C9C1A3"
  )) +
  labs(
    title = "Training Population Heritability (H²)",
    y = "Heritability (H²)",
    x = NULL
  ) +
  theme_minimal(base_size = 14)


# Narrow-sense heritabilities for training population
library(sommer)
library(tidyverse)

# ============================================================
# 1. Filter training population to genotypes with ≥ 2 observations
# ============================================================

train_counts <- train2 %>%
  count(ID, name = "n_obs")

train_filt <- train2 %>%
  inner_join(train_counts %>% filter(n_obs >= 2), by = "ID")

# Subset GRM to matching IDs
keep_ids <- unique(train_filt$ID)
G_train <- G2[keep_ids, keep_ids]

# ============================================================
# 2. Function to compute narrow-sense h2 using sommer
# ============================================================

calc_h2_narrow_train <- function(data, trait, Gmat) {
  df <- data %>% filter(!is.na(.data[[trait]]))
  n_obs <- nrow(df)

  mod <- mmer(
    fixed  = as.formula(paste0(trait, " ~ 1")),
    random = ~ vs(ID, Gu = Gmat),
    data   = df
  )

  vc <- summary(mod)$varcomp

  # Identify additive and residual rows
  row_A <- grep("u:ID", rownames(vc), ignore.case = TRUE)[1]
  row_E <- grep("units", rownames(vc), ignore.case = TRUE)[1]

  sigma_A <- vc[row_A, "VarComp"]
  sigma_E <- vc[row_E, "VarComp"]

  h2 <- sigma_A / (sigma_A + sigma_E)
  SE <- sqrt(2 * (h2^2) / n_obs)

  tibble(
    Trait = trait,
    h2 = h2,
    SE = SE
  )
}

# ============================================================
# 3. Run for all traits
# ============================================================

training_h2_narrow <- map_df(
  c("DON", "INC", "SEV"),
  ~ calc_h2_narrow_train(train_filt, .x, G_train)
) %>%
  mutate(Population = "Training")

training_h2_narrow

ggplot(training_h2_narrow, aes(x = Trait, y = h2, fill = Trait)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = h2 - SE, ymax = h2 + SE),
    width = 0.15,
    linewidth = 0.8
  ) +
  scale_fill_manual(values = c(
    DON = "#6A8F72",
    INC = "#D6B36A",
    SEV = "#C9C1A3"
  )) +
  labs(
    title = "Training Population Narrow-sense Heritability (h²)",
    y = "Heritability (h²)",
    x = NULL
  ) +
  theme_minimal(base_size = 14)
