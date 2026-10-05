read_phenotypes <- function(path, population, environment = NULL) {
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = TRUE, fileEncoding = "UTF-8-BOM")
  aliases <- c(INC = "FHB.incidence.....CO_321.0001149", SEV = "FHB.severity.....CO_321.0001440",
               DON = "FHB.DON.content...ppm.CO_321.0001154", FDK = "FHB.grain.incidence.....CO_321.0001155")
  if (population == "Historical") {
    for (trait in names(aliases)) if (aliases[[trait]] %in% names(x)) names(x)[names(x) == aliases[[trait]]] <- trait
    required_columns(x, c("germplasmName", "FullSampleName", "replicate", "studyName"), path)
    x$genotype <- as.character(x$germplasmName)
    x$REP <- x$replicate
    x$environment <- paste0("Historical::", x$studyName)
  } else {
    required_columns(x, c("ID", "REP"), path)
    x$genotype <- as.character(x$ID)
    if (population == "2025") {
      required_columns(x, c("FullSampleName", "SUB_NURNAME", "EXPT", "YR"), path)
      x <- x[toupper(trimws(x$EXPT)) == "UIUC" & x$YR == 2025, ]
      x$environment <- paste0("2025::", x$SUB_NURNAME)
    } else {
      x$FullSampleName <- x$genotype
      x$environment <- environment
      names(x)[names(x) == "Incidence"] <- "INC"
      names(x)[names(x) == "Severity"] <- "SEV"
      x <- x[x$genotype != "Triticale", ]
    }
  }
  if (anyNA(x$genotype) || any(!nzchar(trimws(x$genotype)))) stop("Missing genotype IDs: ", path)
  for (trait in intersect(config$traits, names(x))) {
    original <- x[[trait]]
    x[[trait]] <- suppressWarnings(as.numeric(original))
    if (any(!is.na(original) & is.na(x[[trait]]))) stop("Non-numeric ", trait, " in ", path)
  }
  x$population <- population
  x
}
long_phenotypes <- function(x) {
  tidyr::pivot_longer(x, cols = dplyr::all_of(intersect(config$traits, names(x))),
                      names_to = "TRAIT", values_to = "y")
}
prepare_analysis <- function() {
  historical <- read_phenotypes(config$historical, "Historical")
  testing <- read_phenotypes(config$testing_2025, "2025")
  locations <- lapply(names(config$locations), function(loc) read_phenotypes(config$locations[[loc]], "2026", loc))
  names(locations) <- names(config$locations)
  required_columns(locations$UNL, "CATEGORY", "UNL validation selection")
  validation_names <- unique(locations$UNL$genotype[locations$UNL$CATEGORY == "Duplicate"])
  if (length(validation_names) != 60L || anyNA(validation_names)) stop("Expected 60 replicated UNL validation genotypes.")

  vcf <- gaston::read.vcf(config$markers, convert.chr = FALSE)
  keep <- is.finite(vcf@snps$maf) & vcf@snps$maf >= config$maf &
    is.finite(vcf@snps$hz) & vcf@snps$hz <= config$heterozygosity
  if (!any(keep)) stop("No markers retained.")
  markers <- gaston::as.matrix(vcf[, keep])
  if (is.null(rownames(markers)) || anyDuplicated(rownames(markers))) stop("Invalid marker sample IDs.")
  # A.mat performs marker-mean imputation; retain that same matrix for PCA.
  relationship <- sommer::A.mat(markers, return.imputed = TRUE)
  grm <- relationship$A
  markers <- relationship$X
  if (any(!is.finite(markers)) || any(!is.finite(grm))) stop("Non-finite marker/GRM values.")
  if (!identical(rownames(grm), colnames(grm)) || !isTRUE(all.equal(grm, t(grm)))) stop("Invalid GRM.")
  historical <- harmonize(historical, grm)
  testing <- harmonize(testing, grm)
  locations <- lapply(locations, harmonize, grm = grm)
  check_mapping(historical, "genotype")
  check_mapping(testing, "genotype")
  selected <- harmonize(data.frame(FullSampleName = validation_names), grm)$FullSampleName
  if (anyDuplicated(selected) || !all(selected %in% rownames(grm))) stop("Some validation genotypes do not match the GRM.")
  for (loc in names(locations)) {
    if (!all(selected %in% locations[[loc]]$FullSampleName)) stop("Missing validation entries at ", loc)
  }

  historical_se <- single_estimates(long_phenotypes(historical))
  testing_se <- single_estimates(long_phenotypes(testing[!testing$SUB_NURNAME %in% config$unreplicated_2025, ]))
  evaluation_se <- single_estimates(long_phenotypes(dplyr::bind_rows(locations)))
  historical_me <- across_estimates(historical_se, "historical across environments")
  evaluation_me <- across_estimates(evaluation_se, "2026 across locations")
  eligible <- function(x) x[!x$genotype %in% validation_names & !x$FullSampleName %in% selected, ]
  expanded_se <- dplyr::bind_rows(eligible(historical_se), eligible(testing_se))
  expanded_me <- across_estimates(expanded_se, "expanded across environments")
  if (length(intersect(selected, expanded_me$FullSampleName))) stop("Expanded validation leakage.")

  objects <- list(GRM = grm, markers = markers, historical_se = historical_se, historical_me = historical_me,
                  testing_se = testing_se, evaluation_se = evaluation_se, evaluation_me = evaluation_me,
                  expanded_me = expanded_me, validation_ids = selected,
                  raw = list(historical = historical, testing = testing, locations = locations))
  for (name in names(objects)) save_result(objects[[name]], paste0(name, ".rds"))
  write_table(data.frame(samples = nrow(markers), markers_input = length(keep), markers_retained = sum(keep),
                         maf_min = config$maf, heterozygosity_max = config$heterozygosity), "marker_summary.csv")
  coverage <- dplyr::summarise(dplyr::group_by(evaluation_se, environment, TRAIT),
                              n_genotypes = dplyr::n_distinct(FullSampleName),
                              n_validation = sum(FullSampleName %in% selected), .groups = "drop")
  write_table(coverage, "phenotype_coverage_2026.csv")
  save_result(preparation_signature(), "preparation_manifest.rds")
}
