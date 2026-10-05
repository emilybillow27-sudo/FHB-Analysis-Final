# Genomic prediction of Fusarium head blight resistance
# Emily Billow
# September 2026

# This script calculates single- and multi-environment BLUEs for the
# historical, 2025, and 2026 populations and evaluates forward genomic
# prediction using GBLUP.

# Run from the project root (FHB Analysis Final).
if (!file.exists("data/FHB_Project_Training_Data.csv")) {
  stop("Set the working directory to the FHB Analysis Final project root.")
}
for (output_dir in c("results/intermediate", "results/tables", "results/figures")) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

library(dplyr)
library(tidyr)
library(purrr)
library(sommer)

traits <- c("INC", "SEV", "DON", "FDK")

# Read in the 2025 testing and historical training datasets

test_raw <- read.csv(
  "data/FHB_Project_Testing_Data_2025.csv"
)

train_raw <- read.csv(
  "data/FHB_Project_Training_Data.csv"
)

test_names <- test_raw %>%
  dplyr::distinct(ID, FullSampleName)

train_names <- train_raw %>%
  dplyr::distinct(germplasmName, FullSampleName)

# Rename the historical trait columns without overwriting train_raw

train_hist <- train_raw %>%
  dplyr::rename(
    DON = `FHB.DON.content...ppm.CO_321.0001154`,
    INC = `FHB.incidence.....CO_321.0001149`,
    SEV = `FHB.severity.....CO_321.0001440`,
    FDK = `FHB.grain.incidence.....CO_321.0001155`
  ) %>%
  dplyr::select(
    studyYear,
    programName,
    studyName,
    studyDescription,
    studyDesign,
    locationName,
    uID,
    FullSampleName,
    germplasmName,
    replicate,
    plotNumber,
    dplyr::all_of(traits)
  )

# Estimate BLUEs for the 2025 testing population
# Analyze each replicated nursery separately

test_blue_list <- list()

for (nursery in unique(test_raw$SUB_NURNAME)) {
  message("Processing 2025 nursery: ", nursery)

  dat_nursery <- test_raw %>%
    dplyr::filter(SUB_NURNAME == nursery)

  if (nursery %in% c("Q Qual AYN", "Topcross", "DH")) {
    warning(
      nursery,
      " is not replicated; retaining it only for raw summaries."
    )

    next
  }

  for (trait in traits) {
    message("Processing trait: ", trait)

    dat_trait <- dat_nursery %>%
      dplyr::select(
        ID,
        FullSampleName,
        REP,
        y = dplyr::all_of(trait)
      ) %>%
      tidyr::drop_na(y)

    if (nrow(dat_trait) == 0) {
      warning(
        "No observations for ",
        trait,
        " in ",
        nursery
      )

      next()
    }

    mod <- sommer::mmes(
      fixed = y ~ ID,
      random = ~REP,
      rcov = ~units,
      data = dat_trait,
      dateWarning = FALSE
    )

    pred <- sommer::predict.mmes(
      mod,
      D = "ID"
    )$pvals %>%
      dplyr::mutate(
        TRAIT = trait,
        SUB_NURNAME = nursery
      ) %>%
      dplyr::left_join(
        test_names,
        by = "ID"
      ) %>%
      dplyr::select(
        ID,
        FullSampleName,
        TRAIT,
        SUB_NURNAME,
        predicted.value,
        std.error
      ) %>%
      dplyr::mutate(
        adjusted = TRUE
      )

    test_blue_list[[length(test_blue_list) + 1]] <- pred
  }
}

blues_se_test <- dplyr::bind_rows(
  test_blue_list
)

saveRDS(
  blues_se_test,
  "results/intermediate/blues_se_test.rds"
)

# Estimate BLUEs for the historical training population
# Process one study and one trait at a time

train_blue_list <- list()

for (study in unique(train_hist$studyName)) {
  message("Processing historical study: ", study)

  dat_study <- train_hist %>%
    dplyr::filter(studyName == study)

  for (trait in traits) {
    message("Processing trait: ", trait)

    dat_trait <- dat_study %>%
      dplyr::select(
        uID,
        FullSampleName,
        germplasmName,
        replicate,
        y = dplyr::all_of(trait)
      ) %>%
      tidyr::drop_na(y)

    if (nrow(dat_trait) == 0) {
      warning(
        "No observations for ",
        trait,
        " in ",
        study
      )

      next()
    }

    if (dplyr::n_distinct(dat_trait$replicate) == 1) {
      warning(
        "There is one replicate for ",
        trait,
        " in ",
        study
      )

      pred <- dat_trait %>%
        dplyr::transmute(
          germplasmName,
          FullSampleName,
          TRAIT = trait,
          studyName = study,
          predicted.value = y,
          std.error = NA_real_,
          adjusted = FALSE
        )

      train_blue_list[[length(train_blue_list) + 1]] <- pred

      next
    }

    mod <- sommer::mmes(
      fixed = y ~ germplasmName,
      random = ~replicate,
      rcov = ~units,
      data = dat_trait,
      dateWarning = FALSE
    )

    pred <- sommer::predict.mmes(
      mod,
      D = "germplasmName"
    )$pvals %>%
      dplyr::mutate(
        TRAIT = trait,
        studyName = study,
        adjusted = TRUE
      ) %>%
      dplyr::left_join(
        train_names,
        by = "germplasmName"
      ) %>%
      dplyr::select(
        germplasmName,
        FullSampleName,
        TRAIT,
        studyName,
        predicted.value,
        std.error,
        adjusted
      )

    train_blue_list[[length(train_blue_list) + 1]] <- pred
  }
}

blues_se_train <- dplyr::bind_rows(
  train_blue_list
)

saveRDS(
  blues_se_train,
  "results/intermediate/blues_se_train.rds"
)

# Estimate precision-weighted historical genotype performance
# across environments

across_train_list <- list()

for (trait in traits) {
  message(
    "Processing across-environment trait: ",
    trait
  )

  dat_trait <- blues_se_train %>%
    dplyr::filter(
      TRAIT == trait,
      !is.na(predicted.value)
    ) %>%
    dplyr::transmute(
      germplasmName,
      FullSampleName,
      studyName,
      y = predicted.value,
      se = std.error
    )

  if (nrow(dat_trait) == 0) {
    warning(
      "No observations for ",
      trait,
      " across environments"
    )

    next()
  }

  if (dplyr::n_distinct(dat_trait$studyName) == 1) {
    pred <- dat_trait %>%
      dplyr::group_by(
        germplasmName,
        FullSampleName
      ) %>%
      dplyr::summarise(
        predicted.value = mean(y),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        TRAIT = trait,
        std.error = NA_real_,
        adjusted = FALSE
      ) %>%
      dplyr::select(
        germplasmName,
        FullSampleName,
        TRAIT,
        predicted.value,
        std.error,
        adjusted
      )

    across_train_list[[length(across_train_list) + 1]] <- pred

    next
  }

  check <- dat_trait %>%
    tidyr::drop_na(se)

  if (nrow(check) == 0) {
    mod <- sommer::mmes(
      fixed = y ~ germplasmName,
      random = ~studyName,
      rcov = ~units,
      data = dat_trait,
      dateWarning = FALSE
    )
  } else if (
    nrow(check) != 0 &
      nrow(check) != nrow(dat_trait)
  ) {
    dat_trait <- dat_trait %>%
      dplyr::mutate(
        se_avg = ifelse(
          is.na(se),
          mean(se, na.rm = TRUE),
          se
        )
      )

    mod <- sommer::mmes(
      fixed = y ~ germplasmName,
      random = ~studyName,
      rcov = ~units,
      W = diag(1 / dat_trait$se_avg^2),
      data = dat_trait,
      dateWarning = FALSE
    )
  } else {
    mod <- sommer::mmes(
      fixed = y ~ germplasmName,
      random = ~studyName,
      rcov = ~units,
      W = diag(1 / dat_trait$se^2),
      data = dat_trait,
      dateWarning = FALSE
    )
  }

  pred <- sommer::predict.mmes(
    mod,
    D = "germplasmName"
  )$pvals %>%
    dplyr::mutate(
      TRAIT = trait,
      adjusted = TRUE
    ) %>%
    dplyr::left_join(
      train_names,
      by = "germplasmName"
    ) %>%
    dplyr::select(
      germplasmName,
      FullSampleName,
      TRAIT,
      predicted.value,
      std.error,
      adjusted
    )

  across_train_list[[length(across_train_list) + 1]] <- pred
}

blues_me_train <- dplyr::bind_rows(
  across_train_list
)

saveRDS(
  blues_me_train,
  "results/intermediate/blues_me_train.rds"
)

# Read and prepare the genome-wide marker data

vcf <- gaston::read.vcf(
  file = "data/fhb_analysis_2026_production_final.vcf.gz",
  convert.chr = FALSE
)

vcf_snps <- vcf@snps %>%
  dplyr::filter(
    maf >= 0.05 &
      hz <= 0.10
  )

vcf <- gaston::as.matrix(
  vcf[
    ,
    vcf@snps$id %in% vcf_snps$id
  ]
)

# Calculate the genomic relationship matrix

grm_inputs <- sommer::A.mat(vcf, return.imputed = TRUE)
GRM <- grm_inputs$A
geno_mat <- as.data.frame(grm_inputs$X, check.names = FALSE)
geno_mat <- tibble::rownames_to_column(geno_mat, "FullSampleName")
saveRDS(geno_mat, "results/intermediate/geno_mat.rds")

saveRDS(
  GRM,
  "results/intermediate/GRM.rds"
)

# Harmonize phenotype identifiers with GRM identifiers

clean_id <- function(x) {
  x %>%
    stringr::str_to_upper() %>%
    stringr::str_squish() %>%
    stringr::str_replace_all(
      "[^A-Z0-9]",
      ""
    )
}

marker_lookup <- tibble::tibble(
  FullSampleName = rownames(GRM)
) %>%
  dplyr::transmute(
    id_key = clean_id(FullSampleName),
    marker_id = FullSampleName
  ) %>%
  dplyr::distinct()

ambiguous_marker_keys <- marker_lookup %>%
  dplyr::count(id_key) %>%
  dplyr::filter(n > 1)

if (nrow(ambiguous_marker_keys) > 0) {
  stop(
    "Some normalized marker IDs are ambiguous; resolve them before matching."
  )
}

harmonize_marker_ids <- function(dat) {
  dat %>%
    dplyr::mutate(
      id_key = clean_id(FullSampleName)
    ) %>%
    dplyr::left_join(
      marker_lookup,
      by = "id_key"
    ) %>%
    dplyr::mutate(
      FullSampleName = dplyr::coalesce(
        marker_id,
        FullSampleName
      )
    ) %>%
    dplyr::select(
      -id_key,
      -marker_id
    )
}

blues_me_train <- harmonize_marker_ids(
  blues_me_train
)

blues_se_test <- harmonize_marker_ids(
  blues_se_test
)

saveRDS(blues_me_train, "results/intermediate/blues_me_train.rds")
saveRDS(blues_se_test, "results/intermediate/blues_se_test.rds")

# Forward prediction of the 2025 testing population using GBLUP

forward_2025_list <- list()

for (trait in traits) {
  message(
    "Predicting 2025 trait: ",
    trait
  )

  dat_trait_train <- blues_me_train %>%
    dplyr::filter(
      TRAIT == trait,
      !is.na(FullSampleName),
      !is.na(predicted.value)
    ) %>%
    dplyr::transmute(
      FullSampleName,
      y = predicted.value,
      y_mask = predicted.value,
      population = "Train"
    ) %>%
    dplyr::distinct()

  dat_trait_test <- blues_se_test %>%
    dplyr::filter(
      TRAIT == trait,
      !is.na(FullSampleName),
      !is.na(predicted.value)
    ) %>%
    dplyr::group_by(
      FullSampleName
    ) %>%
    dplyr::summarise(
      y = mean(
        predicted.value,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      y_mask = NA_real_,
      population = "Test"
    )

  dat_trait_train <- dat_trait_train %>%
    dplyr::filter(
      !FullSampleName %in%
        dat_trait_test$FullSampleName
    )

  dat_trait_complete <- dplyr::bind_rows(
    dat_trait_train,
    dat_trait_test
  ) %>%
    dplyr::filter(
      FullSampleName %in%
        rownames(GRM)
    )

  if (nrow(dat_trait_complete) == 0) {
    warning(
      "No marker-aligned observations for ",
      trait
    )

    next()
  }

  ids_complete <- unique(
    dat_trait_complete$FullSampleName
  )

  dat_GRM <- GRM[
    ids_complete,
    ids_complete,
    drop = FALSE
  ]

  dat_trait_complete$FullSampleName <- factor(
    dat_trait_complete$FullSampleName,
    levels = rownames(dat_GRM)
  )

  fit_2025 <- sommer::mmes(
    fixed = y_mask ~ 1,
    random = ~ sommer::vsm(
      sommer::ism(FullSampleName),
      Gu = dat_GRM
    ),
    rcov = ~units,
    data = dat_trait_complete,
    dateWarning = FALSE
  )

  pred_2025 <- sommer::predict.mmes(
    fit_2025,
    D = "FullSampleName"
  )$pvals %>%
    as.data.frame() %>%
    dplyr::transmute(
      FullSampleName = as.character(
        FullSampleName
      ),
      GEBV = predicted.value
    )

  trait_results <- dat_trait_complete %>%
    dplyr::mutate(
      FullSampleName = as.character(
        FullSampleName
      )
    ) %>%
    dplyr::filter(
      population == "Test"
    ) %>%
    dplyr::select(
      FullSampleName,
      BLUE = y
    ) %>%
    dplyr::left_join(
      pred_2025,
      by = "FullSampleName"
    ) %>%
    dplyr::mutate(
      TRAIT = trait
    )

  forward_2025_list[[trait]] <- trait_results
}

forward_2025 <- dplyr::bind_rows(
  forward_2025_list
)

forward_2025_accuracy <- forward_2025 %>%
  dplyr::group_by(
    TRAIT
  ) %>%
  dplyr::summarise(
    N = sum(
      complete.cases(
        BLUE,
        GEBV
      )
    ),
    PREDICTION_ABILITY = cor(
      BLUE,
      GEBV,
      use = "complete.obs"
    ),
    .groups = "drop"
  )

print(
  forward_2025_accuracy
)

write.csv(
  forward_2025,
  "results/tables/forward_predictions_2025.csv",
  row.names = FALSE
)

write.csv(
  forward_2025_accuracy,
  "results/tables/forward_prediction_2025_accuracy.csv",
  row.names = FALSE
)

# Read the 2026 UNL, UIUC, and SDSU testing populations

unl_2026_raw <- read.csv(
  "data/FHB_Project_Testing_Data_UNL_2026.csv",
  fileEncoding = "UTF-8-BOM",
  check.names = FALSE
)

uiuc_2026_raw <- read.csv(
  "data/FHB_Project_Testing_Data_UIUC_2026.csv",
  fileEncoding = "UTF-8-BOM",
  check.names = FALSE
)

sdsu_2026_raw <- read.csv(
  "data/FHB_Project_Testing_Data_SDSU_2026.csv",
  fileEncoding = "UTF-8-BOM",
  check.names = FALSE
)

# Identify the 60 replicated 2026 validation genotypes

selected_2026 <- unl_2026_raw %>%
  dplyr::filter(
    CATEGORY == "Duplicate"
  ) %>%
  dplyr::distinct(ID) %>%
  dplyr::pull(ID)

stopifnot(
  length(selected_2026) == 60
)

print(
  length(selected_2026)
)

# Standardize the 2026 datasets

standardize_2026 <- function(
    dat,
    location
) {
  if (!"REP" %in% names(dat)) {
    stop(
      "The ",
      location,
      " dataset must contain a REP column."
    )
  }
  
  dat %>%
    dplyr::rename(
      INC = Incidence,
      SEV = Severity
    ) %>%
    dplyr::mutate(
      Location = location,
      PLOT = as.character(PLOT),
      ID = as.character(ID),
      FullSampleName = ID,
      REP = as.factor(REP)
    ) %>%
    dplyr::filter(
      ID != "Triticale"
    )
}

unl_2026 <- standardize_2026(
  unl_2026_raw,
  "UNL"
)

uiuc_2026 <- standardize_2026(
  uiuc_2026_raw,
  "UIUC"
)

sdsu_2026 <- standardize_2026(
  sdsu_2026_raw,
  "SDSU"
)

# Identify traits available at each 2026 location

traits_2026_unl <- traits[
  traits %in% names(unl_2026)
]

traits_2026_uiuc <- traits[
  traits %in% names(uiuc_2026)
]

traits_2026_sdsu <- traits[
  traits %in% names(sdsu_2026)
]

print(
  traits_2026_unl
)

print(
  traits_2026_uiuc
)

print(
  traits_2026_sdsu
)

# Identify all traits currently available in the 2026 dataset

traits_2026 <- unique(
  c(
    traits_2026_unl,
    traits_2026_uiuc,
    traits_2026_sdsu
  )
)

print(
  traits_2026
)

# Calculate single-environment BLUEs for the 2026 populations

calculate_2026_estimates <- function(
    dat,
    location,
    location_traits
) {
  output <- list()
  
  for (trait in location_traits) {
    message(
      "Processing ",
      location,
      " 2026 trait: ",
      trait
    )
    
    dat_trait <- dat %>%
      dplyr::select(
        ID,
        FullSampleName,
        REP,
        y = dplyr::all_of(trait)
      ) %>%
      tidyr::drop_na(y)
    
    if (nrow(dat_trait) == 0) {
      warning(
        "No observations for ",
        trait,
        " at ",
        location
      )
      
      next()
    }
    
    if (
      nrow(dat_trait) ==
      dplyr::n_distinct(dat_trait$ID)
    ) {
      estimates <- dat_trait %>%
        dplyr::transmute(
          ID,
          FullSampleName,
          TRAIT = trait,
          Location = location,
          predicted.value = y,
          std.error = NA_real_,
          adjusted = FALSE
        )
      
      output[[length(output) + 1]] <- estimates
      
      next()
    }
    
    model <- sommer::mmes(
      fixed = y ~ ID,
      random = ~REP,
      rcov = ~units,
      data = dat_trait,
      dateWarning = FALSE
    )
    
    estimates <- sommer::predict.mmes(
      model,
      D = "ID"
    )$pvals %>%
      as.data.frame() %>%
      dplyr::transmute(
        ID = as.character(ID),
        FullSampleName = ID,
        TRAIT = trait,
        Location = location,
        predicted.value = predicted.value,
        std.error = std.error,
        adjusted = TRUE
      )
    
    output[[length(output) + 1]] <- estimates
  }
  
  dplyr::bind_rows(
    output
  )
}

blues_2026_unl <- calculate_2026_estimates(
  unl_2026,
  "UNL",
  traits_2026_unl
)

blues_2026_uiuc <- calculate_2026_estimates(
  uiuc_2026,
  "UIUC",
  traits_2026_uiuc
)

blues_2026_sdsu <- calculate_2026_estimates(
  sdsu_2026,
  "SDSU",
  traits_2026_sdsu
)

blues_se_2026 <- dplyr::bind_rows(
  blues_2026_unl,
  blues_2026_uiuc,
  blues_2026_sdsu
)

# Harmonize 2026 identifiers with GRM identifiers

blues_se_2026 <- harmonize_marker_ids(
  blues_se_2026
)

selected_2026_lookup <- tibble::tibble(
  FullSampleName = selected_2026
)

selected_2026_lookup <- harmonize_marker_ids(
  selected_2026_lookup
)

selected_2026_marker_ids <- selected_2026_lookup %>%
  dplyr::pull(
    FullSampleName
  ) %>%
  unique()

stopifnot(
  length(selected_2026_marker_ids) == 60
)

# Summarize the 2026 single-environment BLUEs

phenotype_summary_2026 <- blues_se_2026 %>%
  dplyr::group_by(
    Location,
    TRAIT
  ) %>%
  dplyr::summarise(
    N = dplyr::n(),
    N_GENOTYPES = dplyr::n_distinct(
      FullSampleName
    ),
    N_VALIDATION_GENOTYPES = dplyr::n_distinct(
      FullSampleName[
        FullSampleName %in%
          selected_2026_marker_ids
      ]
    ),
    .groups = "drop"
  )

print(
  phenotype_summary_2026
)

saveRDS(
  blues_se_2026,
  "results/intermediate/blues_se_2026.rds"
)

# Calculate precision-weighted multi-environment BLUEs for 2026

across_2026_list <- list()

for (trait in traits_2026) {
  message(
    "Processing 2026 across-environment trait: ",
    trait
  )
  
  dat_trait <- blues_se_2026 %>%
    dplyr::filter(
      TRAIT == trait,
      !is.na(predicted.value),
      !is.na(FullSampleName)
    ) %>%
    dplyr::transmute(
      FullSampleName,
      Location,
      y = predicted.value,
      se = std.error
    )
  
  if (nrow(dat_trait) == 0) {
    warning(
      "No 2026 observations for ",
      trait
    )
    
    next()
  }
  
  if (dplyr::n_distinct(dat_trait$Location) == 1) {
    pred <- dat_trait %>%
      dplyr::group_by(
        FullSampleName
      ) %>%
      dplyr::summarise(
        predicted.value = mean(y),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        TRAIT = trait,
        std.error = NA_real_,
        adjusted = FALSE
      ) %>%
      dplyr::select(
        FullSampleName,
        TRAIT,
        predicted.value,
        std.error,
        adjusted
      )
    
    across_2026_list[[trait]] <- pred
    
    next
  }
  
  dat_trait$FullSampleName <- factor(
    dat_trait$FullSampleName
  )
  
  dat_trait$Location <- factor(
    dat_trait$Location
  )
  
  check <- dat_trait %>%
    tidyr::drop_na(se)
  
  if (nrow(check) == 0) {
    mod <- sommer::mmes(
      fixed = y ~ FullSampleName,
      random = ~Location,
      rcov = ~units,
      data = dat_trait,
      dateWarning = FALSE
    )
  } else if (
    nrow(check) != 0 &
    nrow(check) != nrow(dat_trait)
  ) {
    dat_trait <- dat_trait %>%
      dplyr::mutate(
        se_avg = ifelse(
          is.na(se),
          mean(se, na.rm = TRUE),
          se
        )
      )
    
    mod <- sommer::mmes(
      fixed = y ~ FullSampleName,
      random = ~Location,
      rcov = ~units,
      W = diag(1 / dat_trait$se_avg^2),
      data = dat_trait,
      dateWarning = FALSE
    )
  } else {
    mod <- sommer::mmes(
      fixed = y ~ FullSampleName,
      random = ~Location,
      rcov = ~units,
      W = diag(1 / dat_trait$se^2),
      data = dat_trait,
      dateWarning = FALSE
    )
  }
  
  pred <- sommer::predict.mmes(
    mod,
    D = "FullSampleName"
  )$pvals %>%
    as.data.frame() %>%
    dplyr::transmute(
      FullSampleName = as.character(
        FullSampleName
      ),
      TRAIT = trait,
      predicted.value = predicted.value,
      std.error = std.error,
      adjusted = TRUE
    )
  
  across_2026_list[[trait]] <- pred
}

blues_me_2026_all <- dplyr::bind_rows(
  across_2026_list
)

# Restrict the primary forward-validation phenotype set to
# the 60 replicated genotypes

blues_me_2026 <- blues_me_2026_all %>%
  dplyr::filter(
    FullSampleName %in%
      selected_2026_marker_ids
  )

blue_summary_2026 <- blues_me_2026 %>%
  dplyr::group_by(
    TRAIT
  ) %>%
  dplyr::summarise(
    N = dplyr::n(),
    N_GENOTYPES = dplyr::n_distinct(
      FullSampleName
    ),
    .groups = "drop"
  )

print(
  blue_summary_2026
)

saveRDS(
  blues_me_2026_all,
  "results/intermediate/blues_me_2026_all.rds"
)

saveRDS(
  blues_me_2026,
  "results/intermediate/blues_me_2026_validation.rds"
)

# Load the expanded training population calculated in the
# expanded-training analysis

blues_me_expanded <- readRDS(
  "results/intermediate/blues_me_expanded.rds"
)

blues_me_expanded <- harmonize_marker_ids(
  blues_me_expanded
)

# Confirm that the 2026 validation genotypes were excluded
# from the expanded training population

validation_overlap <- intersect(
  selected_2026_marker_ids,
  unique(
    blues_me_expanded$FullSampleName
  )
)

stopifnot(
  length(validation_overlap) == 0
)

# Predict 2026 performance using GBLUP

predict_2026_gblup <- function(
  train_pheno,
  target_pheno,
  trait,
  GRM,
  scenario
) {
  ph_train <- train_pheno %>%
    dplyr::filter(
      TRAIT == trait,
      !is.na(FullSampleName),
      !is.na(predicted.value)
    ) %>%
    dplyr::transmute(
      FullSampleName,
      y = predicted.value,
      y_mask = predicted.value,
      population = "Train"
    ) %>%
    dplyr::distinct()

  ph_target <- target_pheno %>%
    dplyr::filter(
      TRAIT == trait,
      !is.na(FullSampleName),
      !is.na(predicted.value)
    ) %>%
    dplyr::transmute(
      FullSampleName,
      y = predicted.value,
      y_mask = NA_real_,
      population = "Test"
    ) %>%
    dplyr::distinct()

  ph_train <- ph_train %>%
    dplyr::filter(
      !FullSampleName %in%
        ph_target$FullSampleName
    )

  ph_train <- ph_train %>%
    dplyr::filter(
      FullSampleName %in%
        rownames(GRM)
    )

  ph_target <- ph_target %>%
    dplyr::filter(
      FullSampleName %in%
        rownames(GRM)
    )

  if (
    nrow(ph_train) < 3 ||
      nrow(ph_target) == 0
  ) {
    return(NULL)
  }

  dat_complete <- dplyr::bind_rows(
    ph_train,
    ph_target
  )

  ids_complete <- unique(
    dat_complete$FullSampleName
  )

  GRM_trait <- GRM[
    ids_complete,
    ids_complete,
    drop = FALSE
  ]

  dat_complete$FullSampleName <- factor(
    dat_complete$FullSampleName,
    levels = rownames(GRM_trait)
  )

  fit <- sommer::mmes(
    fixed = y_mask ~ 1,
    random = ~ sommer::vsm(
      sommer::ism(FullSampleName),
      Gu = GRM_trait
    ),
    rcov = ~units,
    data = dat_complete,
    dateWarning = FALSE
  )

  pred <- sommer::predict.mmes(
    fit,
    D = "FullSampleName"
  )$pvals %>%
    as.data.frame() %>%
    dplyr::transmute(
      FullSampleName = as.character(
        FullSampleName
      ),
      GEBV = predicted.value
    )

  dat_complete %>%
    dplyr::mutate(
      FullSampleName = as.character(
        FullSampleName
      )
    ) %>%
    dplyr::filter(
      population == "Test"
    ) %>%
    dplyr::select(
      FullSampleName,
      BLUE = y
    ) %>%
    dplyr::left_join(
      pred,
      by = "FullSampleName"
    ) %>%
    dplyr::mutate(
      TRAIT = trait,
      TRAINING_POPULATION = scenario
    )
}

historical_predictions_2026 <- purrr::map_dfr(
  traits_2026,
  ~ predict_2026_gblup(
    blues_me_train,
    blues_me_2026,
    .x,
    GRM,
    "Historical"
  )
)

expanded_predictions_2026 <- purrr::map_dfr(
  traits_2026,
  ~ predict_2026_gblup(
    blues_me_expanded,
    blues_me_2026,
    .x,
    GRM,
    "Expanded"
  )
)

predictions_2026 <- dplyr::bind_rows(
  historical_predictions_2026,
  expanded_predictions_2026
)

# Calculate forward prediction ability

results_2026 <- predictions_2026 %>%
  dplyr::group_by(
    TRAINING_POPULATION,
    TRAIT
  ) %>%
  dplyr::summarise(
    N = sum(
      complete.cases(
        BLUE,
        GEBV
      )
    ),
    PREDICTION_ABILITY = cor(
      BLUE,
      GEBV,
      use = "complete.obs"
    ),
    .groups = "drop"
  )

print(
  results_2026
)

# Compare historical and expanded training populations

results_2026_comparison <- results_2026 %>%
  dplyr::select(
    TRAINING_POPULATION,
    TRAIT,
    N,
    PREDICTION_ABILITY
  ) %>%
  tidyr::pivot_wider(
    names_from = TRAINING_POPULATION,
    values_from = c(
      N,
      PREDICTION_ABILITY
    )
  ) %>%
  dplyr::mutate(
    CHANGE_IN_PREDICTION_ABILITY =
      PREDICTION_ABILITY_Expanded -
        PREDICTION_ABILITY_Historical
  )

print(
  results_2026_comparison
)

# Save the 2026 forward-validation results

write.csv(
  predictions_2026,
  "results/tables/forward_predictions_2026.csv",
  row.names = FALSE
)

write.csv(
  results_2026,
  "results/tables/forward_prediction_2026_accuracy.csv",
  row.names = FALSE
)

write.csv(
  results_2026_comparison,
  "results/tables/forward_prediction_2026_comparison.csv",
  row.names = FALSE
)