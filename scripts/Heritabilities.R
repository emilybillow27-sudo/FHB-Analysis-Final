# ============================================================
# Heritabilities.R — uses GRM already in memory from BLUEs.R
# ============================================================

library(tidyverse)
library(lme4)
library(sommer)

# ------------------------------------------------------------
# 0. Load phenotype data
# ------------------------------------------------------------

fhb   <- read_csv("data/FHB_Project_Testing_Data.csv")
train <- read_csv("data/FHB_Project_Training_Data.csv")

# ------------------------------------------------------------
# 1. Confirm GRM exists
# ------------------------------------------------------------

if (!exists("GRM")) stop("GRM not found. Run BLUEs.R first.")
G2 <- GRM

# ------------------------------------------------------------
# 2. Rename training columns safely
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
# 3. TESTING POPULATION — Broad-sense H²
# ------------------------------------------------------------

calc_h2_test <- function(trait) {
  df <- fhb[!is.na(fhb[[trait]]), ]
  n_rep <- length(unique(df$REP))
  n_obs <- nrow(df)
  
  model <- lmer(as.formula(paste0(trait, " ~ (1|ID) + (1|REP)")), data = df)
  vc <- as.data.frame(VarCorr(model))
  
  Vg <- vc$vcov[vc$grp == "ID"]
  Ve <- vc$vcov[vc$grp == "Residual"]
  
  H2 <- Vg / (Vg + Ve / n_rep)
  SE <- sqrt(2 * H2^2 / n_obs)
  
  data.frame(Trait = trait, Pop = "Testing",
             H2 = H2, h2 = NA,
             SE_H2 = SE, SE_h2 = NA,
             Reps = n_rep)
}

# Only include DI in testing if it exists
testing_traits <- c("DON", "INC", "SEV")
if ("DI" %in% colnames(fhb)) testing_traits <- c(testing_traits, "DI")

testing_h2 <- do.call(rbind, lapply(testing_traits, calc_h2_test))

# ------------------------------------------------------------
# 4. TRAINING POPULATION — Broad-sense H²
# ------------------------------------------------------------

calc_h2_train <- function(trait) {
  df <- train[!is.na(train[[trait]]), ]
  n_rep <- length(unique(df$replicate))
  n_obs <- nrow(df)
  
  if (n_rep >= 2) {
    model <- lmer(as.formula(paste0(trait, " ~ (1|ID) + (1|replicate)")), data = df)
  } else {
    model <- lmer(as.formula(paste0(trait, " ~ (1|ID)")), data = df)
  }
  
  vc <- as.data.frame(VarCorr(model))
  Vg <- vc$vcov[vc$grp == "ID"]
  Ve <- vc$vcov[vc$grp == "Residual"]
  
  H2 <- if (n_rep >= 2) Vg / (Vg + Ve / n_rep) else Vg / (Vg + Ve)
  SE <- sqrt(2 * H2^2 / n_obs)
  
  data.frame(Trait = trait, H2 = H2, SE_H2 = SE, Reps = n_rep)
}

# Training traits ALWAYS include DI
training_traits <- c("DON", "INC", "SEV", "DI")

training_h2 <- do.call(rbind, lapply(training_traits, calc_h2_train))

# ------------------------------------------------------------
# 5. TRAINING POPULATION — Narrow-sense h²
# ------------------------------------------------------------

train$FullSampleName <- as.character(train$FullSampleName)

id_tab <- as.data.frame(table(train$FullSampleName))
colnames(id_tab) <- c("FullSampleName", "n_obs")

ids_keep <- id_tab$FullSampleName[id_tab$n_obs >= 2]
train_filt <- train[train$FullSampleName %in% ids_keep, ]

keep_ids <- unique(train_filt$FullSampleName)
G_train <- G2[keep_ids, keep_ids]

calc_h2_narrow_train <- function(trait) {
  df <- train_filt[!is.na(train_filt[[trait]]), ]
  n_obs <- nrow(df)
  
  mod <- mmer(
    fixed  = as.formula(paste0(trait, " ~ 1")),
    random = ~ vs(FullSampleName, Gu = G_train),
    data   = df
  )
  
  vc <- summary(mod)$varcomp
  sigma_A <- vc[grep("FullSampleName", rownames(vc), ignore.case = TRUE)[1], "VarComp"]
  sigma_E <- vc[grep("units", rownames(vc), ignore.case = TRUE)[1], "VarComp"]
  
  h2 <- sigma_A / (sigma_A + sigma_E)
  SE <- sqrt(2 * h2^2 / n_obs)
  
  data.frame(Trait = trait, h2 = h2, SE_h2 = SE)
}

training_h2_narrow <- do.call(rbind, lapply(training_traits, calc_h2_narrow_train))

# ------------------------------------------------------------
# 6. MERGE TRAINING H² + h² INTO SAME ROWS
# ------------------------------------------------------------

training_all <- merge(training_h2, training_h2_narrow, by = "Trait", all = TRUE)
training_all$Pop <- "Training"

training_all <- training_all[, c("Trait","Pop","H2","h2","SE_H2","SE_h2","Reps")]

# ------------------------------------------------------------
# 7. FINAL TABLE + SAVE
# ------------------------------------------------------------

herit_all <- rbind(testing_h2, training_all)

write.csv(herit_all, "results/tables/heritabilities_all.csv", row.names = FALSE)
