# ============================================================
# Variance Components for Training + Testing Populations
# ============================================================

library(tidyverse)
library(lme4)

# ------------------------------------------------------------
# 0. Load phenotype data
# ------------------------------------------------------------

fhb   <- read_csv("data/FHB_Project_Testing_Data.csv")
train <- read_csv("data/FHB_Project_Training_Data.csv")

# ------------------------------------------------------------
# 1. Rename training columns safely
# ------------------------------------------------------------

don_col <- grep("DON.content",   colnames(train), value = TRUE)[1]
inc_col <- grep("incidence",     colnames(train), value = TRUE)[1]
sev_col <- grep("severity",      colnames(train), value = TRUE)[1]
di_col  <- grep("disease.index", colnames(train), value = TRUE)[1]

names(train)[names(train) == don_col] <- "DON"
names(train)[names(train) == inc_col] <- "INC"
names(train)[names(train) == sev_col] <- "SEV"
names(train)[names(train) == di_col]  <- "DI"
names(train)[names(train) == "uID"]   <- "ID"

train$ID <- as.character(train$ID)

# ------------------------------------------------------------
# 2. Define trait lists
# ------------------------------------------------------------

testing_traits  <- c("DON", "INC", "SEV")
if ("DI" %in% colnames(fhb)) testing_traits <- c(testing_traits, "DI")

training_traits <- c("DON", "INC", "SEV", "DI")

# ------------------------------------------------------------
# Helper: safe extractor for VarCorr
# ------------------------------------------------------------

safe_vc <- function(vc, grp) {
  x <- vc$vcov[vc$grp == grp]
  if (length(x) == 0) return(NA)
  return(x)
}

# ------------------------------------------------------------
# 3. TESTING POPULATION — Variance Components (NO G×E)
# ------------------------------------------------------------

calc_varcomps_test <- function(trait) {
  df <- fhb[!is.na(fhb[[trait]]), ]
  
  model <- lmer(as.formula(paste0(trait, " ~ (1|ID) + (1|REP)")), data = df)
  vc <- as.data.frame(VarCorr(model))
  
  Vg   <- safe_vc(vc, "ID")
  Vrep <- safe_vc(vc, "REP")
  Ve   <- safe_vc(vc, "Residual")
  
  Total <- sum(c(Vg, Vrep, Ve), na.rm = TRUE)
  
  data.frame(
    Trait   = trait,
    Pop     = "Testing",
    Vg      = Vg,
    Vge     = NA,
    Vrep    = Vrep,
    Ve      = Ve,
    Total   = Total,
    Prop_G  = Vg   / Total,
    Prop_GE = NA,
    Prop_REP= Vrep / Total,
    Prop_E  = Ve   / Total
  )
}

testing_varcomps <- do.call(rbind, lapply(testing_traits, calc_varcomps_test))

# ------------------------------------------------------------
# 4. TRAINING POPULATION — Variance Components (WITH SAFE G×E)
# ------------------------------------------------------------

calc_varcomps_train <- function(trait) {
  df <- train[!is.na(train[[trait]]), ]
  env_var <- "locationName"
  
  n_id  <- length(unique(df$ID))
  n_env <- length(unique(df[[env_var]]))
  n_rep <- length(unique(df$replicate))
  
  random_terms <- c()
  if (n_id > 1) random_terms <- c(random_terms, "(1|ID)")
  if (n_rep > 1) random_terms <- c(random_terms, "(1|replicate)")
  if (n_id > 1 && n_env > 1) random_terms <- c(random_terms, paste0("(1|ID:", env_var, ")"))
  
  if (length(random_terms) == 0) {
    return(data.frame(
      Trait = trait, Pop = "Training",
      Vg = NA, Vge = NA, Vrep = NA, Ve = NA,
      Total = NA, Prop_G = NA, Prop_GE = NA, Prop_REP = NA, Prop_E = NA
    ))
  }
  
  formula_str <- paste0(trait, " ~ ", paste(random_terms, collapse = " + "))
  model <- lmer(as.formula(formula_str), data = df)
  vc <- as.data.frame(VarCorr(model))
  
  Vg   <- safe_vc(vc, "ID")
  Vge  <- safe_vc(vc, paste0("ID:", env_var))
  Vrep <- safe_vc(vc, "replicate")
  Ve   <- safe_vc(vc, "Residual")
  
  Total <- sum(c(Vg, Vge, Vrep, Ve), na.rm = TRUE)
  
  data.frame(
    Trait   = trait,
    Pop     = "Training",
    Vg      = Vg,
    Vge     = Vge,
    Vrep    = Vrep,
    Ve      = Ve,
    Total   = Total,
    Prop_G  = Vg   / Total,
    Prop_GE = Vge  / Total,
    Prop_REP= Vrep / Total,
    Prop_E  = Ve   / Total
  )
}

training_varcomps <- do.call(rbind, lapply(training_traits, calc_varcomps_train))

# ------------------------------------------------------------
# 5. MERGE BOTH POPULATIONS (ensure G×E columns are included)
# ------------------------------------------------------------

varcomps_all <- rbind(
  testing_varcomps %>%
    dplyr::mutate(Vge = NA, Prop_GE = NA),  # ensure consistent columns
  training_varcomps
)

# Reorder columns for clarity
varcomps_all <- varcomps_all %>%
  dplyr::select(
    Trait, Pop, Vg, Vge, Vrep, Ve, Total,
    Prop_G, Prop_GE, Prop_REP, Prop_E
  )

# Save full table including G×E
write.csv(varcomps_all, "results/tables/variance_components_all.csv", row.names = FALSE)

# ============================================================
# FIGURE WITH G×E INCLUDED
# ============================================================

varcomps_all <- read_csv("results/tables/variance_components_all.csv")

varcomps_filtered <- varcomps_all %>%
  dplyr::filter(!(Pop == "Testing" & Trait == "DI")) %>%
  dplyr::mutate(Trait = factor(Trait, levels = unique(Trait)))

varcomps_long <- varcomps_filtered %>%
  dplyr::select(Trait, Pop, Prop_G, Prop_GE, Prop_REP, Prop_E) %>%
  tidyr::pivot_longer(
    cols = c(Prop_G, Prop_GE, Prop_REP, Prop_E),
    names_to = "Component",
    values_to = "Proportion"
  )

varcomps_long$Component <- factor(
  varcomps_long$Component,
  levels = c("Prop_G", "Prop_GE", "Prop_REP", "Prop_E"),
  labels = c("Genetic (Vg)", "G×E (Vge)", "Replicate (Vrep)", "Residual (Ve)")
)

p <- ggplot(varcomps_long, aes(x = Trait, y = Proportion, fill = Component)) +
  geom_col(color = "black") +
  facet_wrap(~ Pop, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = c(
    "steelblue",
    "orchid",
    "goldenrod",
    "gray70"
  )) +
  labs(
    title = "Proportion of Variance Components by Trait",
    x = "Trait",
    y = "Proportion of Total Variance",
    fill = "Variance Component"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("results/figures/variance_components.png", p, width = 10, height = 6, dpi = 300)
