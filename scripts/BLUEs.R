# Read in phenotypic data
test <- read.csv("data/FHB_Project_Testing_Data.csv")
train <- read.csv("data/FHB_Project_Training_Data.csv")

test_names <- test |>
  dplyr::distinct(ID, FullSampleName)

train_names <- train |>
  dplyr::distinct(germplasmName, FullSampleName)

blues_se_test <- c()

# Calculate adjusted means (testing)
for (j in unique(test$SUB_NURNAME)) {
  subnurname <- test |>
    dplyr::filter(SUB_NURNAME == j)
  if (j %in% c("Q Qual AYN", "Topcross", "DH")) {
    warning(j, " is not replicated. Moving to next nursery.")
    next()
  }

  for (i in c("INC", "SEV", "DON")) {
    print(i)
    temp1 <- subnurname |>
      dplyr::select(ID, FullSampleName, REP, dplyr::all_of(i))
    colnames(temp1)[4] <- "y"
    temp1 <- temp1 |>
      tidyr::drop_na(y)
    temp2 <- sommer::mmes(
      fixed = y ~ ID,
      random = ~REP,
      rcov = ~units,
      data = temp1,
      dateWarning = FALSE
    )
    temp3 <- sommer::predict.mmes(temp2, D = "ID")$pvals
    temp3$TRAIT <- i
    temp3$SUB_NURNAME <- j
    temp3 <- temp3 |>
      dplyr::left_join(test_names, by = "ID") |>
      dplyr::select(ID, FullSampleName, TRAIT, SUB_NURNAME, predicted.value, std.error)
    blues_se_test <- rbind(blues_se_test, temp3)
    rm(temp1, temp2, temp3)
  }
  rm(subnurname)
}

check_test <- blues_se_test |>
  dplyr::distinct(SUB_NURNAME, ID) |>
  dplyr::group_by(ID) |>
  dplyr::count()


# Calculate adjusted means (training)
train <- train |>
  dplyr::rename(
    DON = `FHB.DON.content...ppm.CO_321.0001154`,
    INC = `FHB.incidence.....CO_321.0001149`,
    SEV = `FHB.severity.....CO_321.0001440`
  ) |>
  dplyr::select(
    studyYear, programName, studyName, studyDescription, studyDesign, locationName,
    uID, FullSampleName, germplasmName, replicate, plotNumber,
    INC, SEV, DON
  )

blues_se_train <- c()

for (j in unique(train$studyName)) {
  # Filter for study name
  studyName <- train |>
    dplyr::filter(studyName == j)

  for (i in c("INC", "SEV", "DON")) {
    # Print trait name
    print(i)

    # Select relevant columns
    temp1 <- studyName |>
      dplyr::select(uID, FullSampleName, germplasmName, replicate, dplyr::all_of(i))

    # Rename column 5
    colnames(temp1)[5] <- "y"

    # Drop missing values of y
    temp1 <- temp1 |>
      tidyr::drop_na(y)

    # Check if temp1 has data, if not skip
    if (nrow(temp1) == 0) {
      # Throw warning
      warning("No observations for ", i, " in ", j)

      # Remove temp1
      rm(temp1)

      # Skip to next level of i
      next()
    }

    # Check if there is one replicate
    if (length(unique(temp1$replicate)) == 1) {
      # Warning
      warning("There is one replicate for ", i, " in ", j)

      # Manipulate temp1 to resemble blues_se_train
      temp1 <- temp1 |>
        # Add columns to temp1 to look like temp3
        dplyr::mutate(
          TRAIT = i,
          studyName = j,
          predicted.value = y,
          std.error = NA,
          adjusted = FALSE
        ) |>
        # Selecting only relevant columns
        dplyr::select(germplasmName, FullSampleName, TRAIT, studyName, predicted.value, std.error, adjusted)

      # Bind in BLUEs
      blues_se_train <- rbind(blues_se_train, temp1)

      # Remove temporary object
      rm(temp1)

      # Skip to next level of i
      next()
    }

    # Run mixed model
    temp2 <- sommer::mmes(
      fixed = y ~ germplasmName,
      random = ~replicate,
      rcov = ~units,
      data = temp1,
      dateWarning = FALSE
    )

    # Calculate BLUEs
    temp3 <- sommer::predict.mmes(temp2, D = "germplasmName")$pvals

    # Assign trait in dataframe
    temp3$TRAIT <- i

    # Assign program name
    temp3$studyName <- j

    # Combine by trait name
    temp3 <- temp3 |>
      # Left join on germplasmName
      dplyr::left_join(train_names, by = "germplasmName") |>
      # Selecting only relevant columns
      dplyr::select(germplasmName, FullSampleName, TRAIT, studyName, predicted.value, std.error) |>
      # Add column to indicate if adjusted or not
      dplyr::mutate(adjusted = TRUE)

    # Bind in BLUEs
    blues_se_train <- rbind(blues_se_train, temp3)

    # Remove temporary objects
    rm(temp1, temp2, temp3)
  }
  # Remove studyname object
  rm(studyName)
}

check_train <- blues_se_train |>
  # Exclude replicated names within a location
  dplyr::distinct(studyName, germplasmName) |>
  # Group by germplasmName
  dplyr::group_by(germplasmName) |>
  # Count number of times the germplasmName appears
  dplyr::count()

# Across-environment means for training dataset
train_me <- blues_se_train

blues_me_train <- c()

for (trait in c("INC", "SEV", "DON")) {
  # Print trait name
  message("Processing trait: ", trait)

  # Select relevant columns
  temp <- train_me |>
    dplyr::filter(TRAIT == trait) |>
    dplyr::select(
      germplasmName, FullSampleName,
      studyName,
      y = predicted.value
    ) |>
    tidyr::drop_na(y)

  # Check if data exists
  if (nrow(temp) == 0) {
    warning("No observations for ", trait)
    next()
  }

  # Check replication across environments
  if (length(unique(temp$studyName)) == 1) {
    warning("Only one replicate across environments for ", trait, ". Returning raw means.")

    # Compute raw means
    temp_out <- temp |>
      dplyr::group_by(germplasmName) |>
      dplyr::summarise(predicted.value = mean(y), .groups = "drop") |>
      dplyr::mutate(
        TRAIT = trait,
        std.error = NA,
        adjusted = FALSE
      ) |>
      dplyr::left_join(train_names, by = "germplasmName") |>
      # Match column order
      dplyr::select(
        germplasmName, FullSampleName,
        TRAIT, predicted.value, std.error, adjusted
      )

    # Bind results
    blues_me_train <- rbind(blues_me_train, temp_out)
    next()
  }

  # Fit across-environment model
  mod <- sommer::mmes(
    fixed = y ~ germplasmName + studyName,
    random = ~studyName,
    rcov = ~units,
    data = temp,
    dateWarning = FALSE
  )

  # Predict BLUEs
  pred <- sommer::predict.mmes(mod, D = "germplasmName")$pvals |>
    dplyr::mutate(
      TRAIT = trait,
      adjusted = TRUE
    ) |>
    # Select final columns
    dplyr::left_join(train_names, by = "germplasmName") |>
    dplyr::select(
      germplasmName, FullSampleName,
      TRAIT, predicted.value, std.error, adjusted
    )

  # Bind results
  blues_me_train <- rbind(blues_me_train, pred)
}

# Count environments per genotype
check_me_train <- blues_me_train |>
  dplyr::group_by(germplasmName) |>
  dplyr::count()

# Read in genotype VCF
library(VariantAnnotation)
library(tidyverse)

vcf_file <- "data/fhb_analysis_2026_production_final.vcf.gz"

# Read VCF
vcf <- readVcf(vcf_file, genome = "unknown")

# Extract genotype (GT) matrix
gt <- geno(vcf)$GT # character matrix like "0/0", "0/1", "1/1"

# Convert GT to numeric dosage (0,1,2)
convert_gt <- function(x) {
  x <- gsub("\\|", "/", x) # standardize separators
  case_when(
    x == "0/0" ~ 0,
    x == "0/1" ~ 1,
    x == "1/0" ~ 1,
    x == "1/1" ~ 2,
    TRUE ~ NA_real_
  )
}

geno_mat <- apply(gt, 2, convert_gt) |>
  t() |>
  as.data.frame()

# Add marker names
colnames(geno_mat) <- rownames(gt)

# Add genotype/sample names
geno_mat$genotype <- colnames(gt)

# Reorder so genotype column is first
geno_mat <- geno_mat |> dplyr::relocate(genotype)

# Check dimensions
dim(geno_mat)
head(geno_mat[, 1:6])

# Build Genomic Relationship Matrix (GRM)
library(tidyverse)

# geno_mat must already exist from previous step
# geno_mat: rows = genotypes, columns = markers, values = 0/1/2

# Remove genotype column for matrix operations
geno_numeric <- geno_mat |>
  dplyr::select(-genotype) |>
  as.matrix()

# Impute missing values with marker means
marker_means <- colMeans(geno_numeric, na.rm = TRUE)

for (m in seq_len(ncol(geno_numeric))) {
  geno_numeric[is.na(geno_numeric[, m]), m] <- marker_means[m]
}

# Center genotype matrix (VanRaden Method 1)
p <- marker_means / 2
Z <- sweep(geno_numeric, 2, 2 * p) # subtract 2p

# Compute denominator: 2 * Σ p(1-p)
denom <- 2 * sum(p * (1 - p))

# GRM = ZZ' / denom
GRM <- (Z %*% t(Z)) / denom

# Add row/column names
rownames(GRM) <- geno_mat$genotype
colnames(GRM) <- geno_mat$genotype

# Quick checks
dim(GRM)
GRM[1:5, 1:5]

# Harmonize genotype IDs

# Rename genotype column for clarity
geno_mat <- geno_mat |>
  dplyr::rename(FullSampleName = genotype)

# Check overlap with phenotype IDs
intersect_ids <- intersect(blues_me_train$FullSampleName, geno_mat$FullSampleName)

length(intersect_ids)

# Filter phenotype BLUEs
pheno_use <- blues_me_train |>
  dplyr::filter(FullSampleName %in% intersect_ids)

# Filter genotype matrix
geno_use <- geno_mat |>
  dplyr::filter(FullSampleName %in% intersect_ids)

# Merge phenotype and genotype data
merged_data <- pheno_use |>
  dplyr::left_join(geno_use, by = "FullSampleName")

# Inspect
dim(merged_data)
head(merged_data[, 1:10])

# Align GRM to phenotype order
common_ids <- merged_data$FullSampleName

GRM_aligned <- GRM[common_ids, common_ids]

# Check
dim(GRM_aligned)
GRM_aligned[1:5, 1:5]


# Load required packages
library(dplyr)
library(tidyr)
library(sommer)
library(rrBLUP)

# Rename phenotype column to y for consistency
merged_data <- merged_data %>%
  dplyr::rename(y = predicted.value)

# Define traits to model
traits <- c("INC", "SEV", "DON")

# Initialize list to store rrBLUP GBLUP results
rr_results <- list()

for (t in traits) {
  message("Running rrBLUP GBLUP for trait: ", t)

  # Subset phenotype data for this trait
  ph <- merged_data %>%
    dplyr::filter(TRAIT == t) %>%
    dplyr::select(FullSampleName, y) %>%
    dplyr::filter(!is.na(y))

  # Identify genotypes present in both phenotype data and GRM
  ids_t <- rownames(GRM)[rownames(GRM) %in% ph$FullSampleName]

  # Reorder phenotype data to match GRM order
  ph_t <- ph %>%
    dplyr::filter(FullSampleName %in% ids_t) %>%
    dplyr::arrange(match(FullSampleName, ids_t))

  # Subset GRM to the same genotypes in the same order
  GRM_t <- GRM[ids_t, ids_t]

  # Fit GBLUP model using rrBLUP
  fit <- mixed.solve(
    y = ph_t$y,
    K = GRM_t,
    X = matrix(1, nrow = nrow(ph_t), ncol = 1)
  )

  # Extract GEBVs from the fitted model
  gebv_vec <- as.vector(fit$u)
  names(gebv_vec) <- rownames(GRM_t)

  # Store GEBVs for this trait
  rr_results[[t]] <- data.frame(
    FullSampleName = names(gebv_vec),
    GEBV = gebv_vec,
    TRAIT = t
  )
}

# Combine GEBVs across traits
rr_all <- dplyr::bind_rows(rr_results)


library(rrBLUP)
library(dplyr)
library(tidyr)

# Prepare genotype matrices for training and testing
geno_train <- geno_mat %>%
  dplyr::filter(FullSampleName %in% blues_me_train$FullSampleName)

geno_test <- geno_mat %>%
  dplyr::filter(FullSampleName %in% blues_se_test$FullSampleName)

# Convert genotype data to numeric matrices
M_train <- geno_train %>%
  dplyr::select(-FullSampleName) %>%
  as.matrix()
M_test <- geno_test %>%
  dplyr::select(-FullSampleName) %>%
  as.matrix()

# Initialize list to store forward prediction results
traits <- c("INC", "SEV", "DON")
forward_results <- list()

for (t in traits) {
  message("Forward predicting for trait: ", t)

  # Subset training phenotype data for this trait
  ph_train <- merged_data %>%
    dplyr::filter(TRAIT == t) %>%
    dplyr::select(FullSampleName, y) %>%
    dplyr::filter(!is.na(y))

  # Align training genotype matrix to phenotype order
  ids_train <- ph_train$FullSampleName
  M_train_t <- M_train[match(ids_train, geno_train$FullSampleName), ]

  # Fit rrBLUP marker model
  fit <- mixed.solve(
    y = ph_train$y,
    Z = M_train_t,
    X = matrix(1, nrow = nrow(M_train_t), ncol = 1)
  )

  # Extract marker effects
  marker_effects <- fit$u
  names(marker_effects) <- colnames(M_train_t)

  # Align testing genotype matrix to marker order
  M_test_t <- M_test[, names(marker_effects), drop = FALSE]

  # Compute forward-predicted GEBVs for testing lines
  gebv_test <- as.vector(M_test_t %*% marker_effects)
  names(gebv_test) <- geno_test$FullSampleName

  # Store forward prediction results
  forward_results[[t]] <- data.frame(
    FullSampleName = names(gebv_test),
    GEBV = gebv_test,
    TRAIT = t
  )
}

# Combine forward prediction results across traits
forward_all <- bind_rows(forward_results)


library(dplyr)

# Prepare testing BLUEs by renaming predicted.value to y
test_blues <- blues_se_test %>%
  dplyr::rename(y = predicted.value) %>%
  dplyr::select(FullSampleName, TRAIT, y) %>%
  dplyr::filter(!is.na(y))

# Merge forward-predicted GEBVs with testing BLUEs
forward_eval <- forward_all %>%
  dplyr::left_join(test_blues, by = c("FullSampleName", "TRAIT"))

# Compute prediction accuracy for each trait
prediction_accuracy <- forward_eval %>%
  dplyr::group_by(TRAIT) %>%
  dplyr::summarise(
    n = sum(!is.na(GEBV) & !is.na(y)),
    accuracy = cor(GEBV, y, use = "complete.obs")
  )

prediction_accuracy