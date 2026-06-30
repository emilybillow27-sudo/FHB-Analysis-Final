###############################################
# Load phenotypic data for testing and training
###############################################

test  <- read.csv("data/FHB_Project_Testing_Data.csv")
train <- read.csv("data/FHB_Project_Training_Data.csv")

# Extract unique sample identifiers for merging BLUEs
test_names <- test |>
  dplyr::distinct(ID, FullSampleName)

train_names <- train |>
  dplyr::distinct(germplasmName, FullSampleName)

###############################################
# BLUEs for original testing dataset (INC, SEV, DON)
###############################################

blues_se_test <- c()

for (nur in unique(test$SUB_NURNAME)) {
  
  sub <- test |> dplyr::filter(SUB_NURNAME == nur)
  
  # Skip non‑replicated nurseries
  if (nur %in% c("Q Qual AYN", "Topcross", "DH")) {
    warning(nur, " is not replicated. Skipping.")
    next()
  }
  
  for (trait in c("INC", "SEV", "DON")) {
    
    temp <- sub |>
      dplyr::select(ID, FullSampleName, REP, dplyr::all_of(trait))
    colnames(temp)[4] <- "y"
    
    temp <- temp |> tidyr::drop_na(y)
    
    mod <- sommer::mmes(
      fixed = y ~ ID,
      random = ~ REP,
      rcov   = ~ units,
      data   = temp,
      dateWarning = FALSE
    )
    
    pred <- sommer::predict.mmes(mod, D = "ID")$pvals
    
    pred$TRAIT <- trait
    pred$SUB_NURNAME <- nur
    
    pred <- pred |>
      dplyr::left_join(test_names, by = "ID") |>
      dplyr::select(ID, FullSampleName, TRAIT, SUB_NURNAME,
                    predicted.value, std.error)
    
    blues_se_test <- rbind(blues_se_test, pred)
  }
}

###############################################
# BLUEs for training dataset (single environments)
###############################################

train <- train |>
  dplyr::rename(
    DON = `FHB.DON.content...ppm.CO_321.0001154`,
    INC = `FHB.incidence.....CO_321.0001149`,
    SEV = `FHB.severity.....CO_321.0001440`
  ) |>
  dplyr::select(
    studyYear, programName, studyName, locationName,
    uID, FullSampleName, germplasmName, replicate,
    plotNumber, INC, SEV, DON
  )

blues_se_train <- c()

for (study in unique(train$studyName)) {
  
  sub <- train |> dplyr::filter(studyName == study)
  
  for (trait in c("INC", "SEV", "DON")) {
    
    temp <- sub |>
      dplyr::select(uID, FullSampleName, germplasmName, replicate, dplyr::all_of(trait))
    colnames(temp)[5] <- "y"
    
    temp <- temp |> tidyr::drop_na(y)
    
    if (nrow(temp) == 0) {
      warning("No observations for ", trait, " in ", study)
      next()
    }
    
    if (length(unique(temp$replicate)) == 1) {
      
      warning("One replicate for ", trait, " in ", study)
      
      temp <- temp |>
        dplyr::mutate(
          TRAIT = trait,
          studyName = study,
          predicted.value = y,
          std.error = NA,
          adjusted = FALSE
        ) |>
        dplyr::select(germplasmName, FullSampleName, TRAIT, studyName,
                      predicted.value, std.error, adjusted)
      
      blues_se_train <- rbind(blues_se_train, temp)
      next()
    }
    
    mod <- sommer::mmes(
      fixed = y ~ germplasmName,
      random = ~ replicate,
      rcov   = ~ units,
      data   = temp,
      dateWarning = FALSE
    )
    
    pred <- sommer::predict.mmes(mod, D = "germplasmName")$pvals
    
    pred$TRAIT <- trait
    pred$studyName <- study
    
    pred <- pred |>
      dplyr::left_join(train_names, by = "germplasmName") |>
      dplyr::select(germplasmName, FullSampleName, TRAIT, studyName,
                    predicted.value, std.error) |>
      dplyr::mutate(adjusted = TRUE)
    
    blues_se_train <- rbind(blues_se_train, pred)
  }
}

###############################################
# Across‑environment BLUEs for training dataset
###############################################

train_me <- blues_se_train
blues_me_train <- c()

for (trait in c("INC", "SEV", "DON")) {
  
  temp <- train_me |>
    dplyr::filter(TRAIT == trait) |>
    dplyr::select(germplasmName, FullSampleName, studyName, y = predicted.value) |>
    tidyr::drop_na(y)
  
  if (nrow(temp) == 0) next()
  
  if (length(unique(temp$studyName)) == 1) {
    
    temp_out <- temp |>
      dplyr::group_by(germplasmName) |>
      dplyr::summarise(predicted.value = mean(y), .groups = "drop") |>
      dplyr::mutate(
        TRAIT = trait,
        std.error = NA,
        adjusted = FALSE
      ) |>
      dplyr::left_join(train_names, by = "germplasmName") |>
      dplyr::select(germplasmName, FullSampleName, TRAIT, predicted.value, std.error, adjusted)
    
    blues_me_train <- rbind(blues_me_train, temp_out)
    next()
  }
  
  mod <- sommer::mmes(
    fixed = y ~ germplasmName + studyName,
    random = ~ studyName,
    rcov   = ~ units,
    data   = temp,
    dateWarning = FALSE
  )
  
  pred <- sommer::predict.mmes(mod, D = "germplasmName")$pvals |>
    dplyr::mutate(TRAIT = trait, adjusted = TRUE) |>
    dplyr::left_join(train_names, by = "germplasmName") |>
    dplyr::select(germplasmName, FullSampleName, TRAIT, predicted.value, std.error, adjusted)
  
  blues_me_train <- rbind(blues_me_train, pred)
}

###############################################
# Genotype processing + GRM construction
###############################################

library(VariantAnnotation)
library(tidyverse)

vcf <- readVcf("data/fhb_analysis_2026_production_final.vcf.gz", genome = "unknown")
gt  <- geno(vcf)$GT

convert_gt <- function(x) {
  x <- gsub("\\|", "/", x)
  dplyr::case_when(
    x == "0/0" ~ 0,
    x == "0/1" ~ 1,
    x == "1/0" ~ 1,
    x == "1/1" ~ 2,
    TRUE ~ NA_real_
  )
}

geno_mat <- apply(gt, 2, convert_gt) |> t() |> as.data.frame()
colnames(geno_mat) <- rownames(gt)
geno_mat$FullSampleName <- colnames(gt)

geno_numeric <- geno_mat |> dplyr::select(-FullSampleName) |> as.matrix()

marker_means <- colMeans(geno_numeric, na.rm = TRUE)

for (m in seq_len(ncol(geno_numeric))) {
  geno_numeric[is.na(geno_numeric[, m]), m] <- marker_means[m]
}

p <- marker_means / 2
Z <- sweep(geno_numeric, 2, 2 * p)
denom <- 2 * sum(p * (1 - p))

GRM <- (Z %*% t(Z)) / denom
rownames(GRM) <- geno_mat$FullSampleName
colnames(GRM) <- geno_mat$FullSampleName

###############################################
# Merge phenotypes + genotypes
###############################################

intersect_ids <- intersect(blues_me_train$FullSampleName, geno_mat$FullSampleName)

pheno_use <- blues_me_train |> dplyr::filter(FullSampleName %in% intersect_ids)
geno_use  <- geno_mat       |> dplyr::filter(FullSampleName %in% intersect_ids)

merged_data <- pheno_use |> dplyr::left_join(geno_use, by = "FullSampleName")
merged_data <- merged_data |> dplyr::rename(y = predicted.value)

GRM_aligned <- GRM[merged_data$FullSampleName, merged_data$FullSampleName]

###############################################
# rrBLUP GBLUP (training set)
###############################################

library(rrBLUP)

traits <- c("INC", "SEV", "DON")
rr_results <- list()

for (t in traits) {
  
  ph <- merged_data |> dplyr::filter(TRAIT == t) |> dplyr::select(FullSampleName, y)
  
  ids <- ph$FullSampleName
  GRM_t <- GRM[ids, ids]
  
  fit <- mixed.solve(
    y = ph$y,
    K = GRM_t,
    X = matrix(1, nrow = nrow(GRM_t), ncol = 1)
  )
  
  gebv <- fit$u
  names(gebv) <- ids
  
  rr_results[[t]] <- data.frame(
    FullSampleName = ids,
    GEBV = gebv,
    TRAIT = t
  )
}

rr_all <- dplyr::bind_rows(rr_results)

###############################################
# Forward prediction for original testing dataset
###############################################

geno_train <- geno_mat |> dplyr::filter(FullSampleName %in% merged_data$FullSampleName)
geno_test  <- geno_mat |> dplyr::filter(FullSampleName %in% blues_se_test$FullSampleName)

M_train <- geno_train |> dplyr::select(-FullSampleName) |> as.matrix()
M_test  <- geno_test  |> dplyr::select(-FullSampleName) |> as.matrix()

forward_results <- list()

for (t in traits) {
  
  ph_train <- merged_data |> dplyr::filter(TRAIT == t)
  
  ids_train <- ph_train$FullSampleName
  M_train_t <- M_train[match(ids_train, geno_train$FullSampleName), ]
  
  fit <- mixed.solve(
    y = ph_train$y,
    Z = M_train_t,
    X = matrix(1, nrow = nrow(M_train_t), ncol = 1)
  )
  
  marker_effects <- fit$u
  names(marker_effects) <- colnames(M_train_t)
  
  M_test_t <- M_test[, names(marker_effects), drop = FALSE]
  
  gebv <- as.vector(M_test_t %*% marker_effects)
  names(gebv) <- geno_test$FullSampleName
  
  forward_results[[t]] <- data.frame(
    FullSampleName = names(gebv),
    GEBV = gebv,
    TRAIT = t
  )
}

forward_all <- dplyr::bind_rows(forward_results)

test_blues <- blues_se_test |>
  dplyr::rename(y = predicted.value) |>
  dplyr::select(FullSampleName, TRAIT, y)

forward_eval <- forward_all |>
  dplyr::left_join(test_blues, by = c("FullSampleName", "TRAIT"))

prediction_accuracy <- forward_eval |>
  dplyr::group_by(TRAIT) |>
  dplyr::summarise(
    n = sum(!is.na(GEBV) & !is.na(y)),
    accuracy = cor(GEBV, y, use = "complete.obs")
  )

###############################################
# Cross‑validation (training set)
###############################################

K <- 5
cv_results <- list()

for (t in traits) {
  
  ph <- merged_data |> dplyr::filter(TRAIT == t)
  
  ids <- ph$FullSampleName
  GRM_t <- GRM[ids, ids]
  
  set.seed(123)
  folds <- sample(rep(1:K, length.out = length(ids)))
  
  pred <- rep(NA, length(ids))
  
  for (k in 1:K) {
    
    test_idx <- which(folds == k)
    train_idx <- setdiff(seq_len(length(ids)), test_idx)
    
    y_train <- ph$y
    y_train[test_idx] <- NA
    
    fit <- mixed.solve(
      y = y_train,
      K = GRM_t,
      X = matrix(1, nrow = length(ids), ncol = 1)
    )
    
    pred[test_idx] <- fit$u[test_idx]
  }
  
  cv_results[[t]] <- data.frame(
    TRAIT = t,
    accuracy = cor(pred, ph$y, use = "complete.obs"),
    n = length(ids)
  )
}

cv_summary <- dplyr::bind_rows(cv_results)

###############################################
# 2026 testing dataset: BLUEs + forward prediction
###############################################

test_2026 <- read.csv("data/FHB_Project_Testing_Data_2026.csv")

# Rename columns to match pipeline
test_2026 <- test_2026 |>
  dplyr::rename(
    FullSampleName = Variety,
    REP = PLOT,
    INC = Incidence,
    SEV = Severity
  )

blues_se_test_2026 <- c()

for (nur in unique(test_2026$NURNAME)) {
  
  sub <- test_2026 |> dplyr::filter(NURNAME == nur)
  
  if (nur %in% c("Q Qual AYN", "Topcross", "DH")) next()
  
  for (trait in c("INC", "SEV")) {
    
    temp <- sub |> dplyr::select(FullSampleName, REP, dplyr::all_of(trait))
    colnames(temp)[3] <- "y"
    
    temp <- temp |> tidyr::drop_na(y)
    
    mod <- sommer::mmes(
      fixed = y ~ FullSampleName,
      random = ~ REP,
      rcov   = ~ units,
      data   = temp,
      dateWarning = FALSE
    )
    
    pred <- sommer::predict.mmes(mod, D = "FullSampleName")$pvals
    
    pred$TRAIT <- trait
    pred$NURNAME <- nur
    
    pred <- pred |> dplyr::select(FullSampleName, TRAIT, NURNAME,
                                  predicted.value, std.error)
    
    blues_se_test_2026 <- rbind(blues_se_test_2026, pred)
  }
}

test_blues_2026 <- blues_se_test_2026 |>
  dplyr::rename(y = predicted.value)

geno_test_2026 <- geno_mat |>
  dplyr::filter(FullSampleName %in% test_blues_2026$FullSampleName)

M_test_2026 <- geno_test_2026 |> dplyr::select(-FullSampleName) |> as.matrix()

###############################################
# Marker effects from training data
###############################################

marker_effects_list <- list()

for (t in c("INC", "SEV", "DON")) {
  
  ph_train <- merged_data |> dplyr::filter(TRAIT == t)
  
  ids_train <- ph_train$FullSampleName
  M_train_t <- M_train[match(ids_train, geno_train$FullSampleName), ]
  
  fit <- mixed.solve(
    y = ph_train$y,
    Z = M_train_t,
    X = matrix(1, nrow = nrow(M_train_t), ncol = 1)
  )
  
  marker_effects_list[[t]] <- fit$u
}

###############################################
# Forward prediction for 2026 dataset
###############################################

forward_2026_results <- list()

for (t in c("INC", "SEV")) {
  
  marker_effects <- marker_effects_list[[t]]
  
  M_test_t <- M_test_2026[, names(marker_effects), drop = FALSE]
  
  gebv <- as.vector(M_test_t %*% marker_effects)
  names(gebv) <- geno_test_2026$FullSampleName
  
  forward_2026_results[[t]] <- data.frame(
    FullSampleName = names(gebv),
    GEBV = gebv,
    TRAIT = t
  )
}

forward_2026_all <- dplyr::bind_rows(forward_2026_results)

forward_eval_2026 <- forward_2026_all |>
  dplyr::left_join(test_blues_2026, by = c("FullSampleName", "TRAIT"))

prediction_accuracy_2026 <- forward_eval_2026 |>
  dplyr::group_by(TRAIT) |>
  dplyr::summarise(
    n = sum(!is.na(GEBV) & !is.na(y)),
    accuracy = cor(GEBV, y, use = "complete.obs")
  )
