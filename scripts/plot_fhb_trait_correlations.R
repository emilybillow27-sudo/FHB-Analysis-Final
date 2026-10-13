# Standalone descriptive correlations; no models are refitted.
# UIUC observations are averaged within ID and sub-nursery.
# Only replicated Elite and WSS AYN sub-nurseries are included; sub-nurseries remain separate.
# These descriptive correlations do not measure genomic predictive ability.
library(dplyr)
library(tidyr)
project_dir <- getOption('fhb.project_dir', '/Users/emilybillow/Desktop/FHB Analysis Final')
out_dir <- getOption('fhb.table_dir', file.path(project_dir, 'results', 'tables', 'trait_correlations'))
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
traits <- c('INC','SEV','DON','FDK')
labels <- c(INC='INC', SEV='SEV', DON='DON', FDK='FDK')
prepare <- function(file, id, uiuc=FALSE) {
  candidates <- file.path(project_dir, c('results/intermediate','results'), file)
  available <- candidates[file.exists(candidates)]
  if (!length(available)) stop('Missing input: ', file)
  if (length(available)>1) stop('Two copies of ', file, ' found. Specify the intended input before proceeding.')
  message('Reading: ', available)
  d <- readRDS(available)
  stopifnot(all(c(id,'TRAIT','predicted.value','adjusted') %in% names(d)))
  d <- d |> filter(TRAIT %in% traits, adjusted == TRUE)
  if (uiuc) {
    stopifnot('SUB_NURNAME' %in% names(d))
    if (any(is.na(d$SUB_NURNAME)) || !all(d$SUB_NURNAME %in% c('Elite','WSS AYN')))
      stop('Check UIUC BLUE sub-nursery membership.')
    keys <- c(id,'SUB_NURNAME','TRAIT')
  } else keys <- c(id,'TRAIT')
  if (anyDuplicated(d[keys])) stop('Duplicate estimates in ', file)
  if (any(is.na(d[[id]]) | trimws(as.character(d[[id]]))=='')) stop('Missing entry identifiers')
  if (!is.numeric(d$predicted.value)) stop('Estimates must be numeric')
  if (!all(traits %in% d$TRAIT)) stop('An expected trait is missing')
  w <- d |> group_by(across(all_of(id)),TRAIT) |>
    summarise(value=if(all(!is.finite(predicted.value))) NA_real_ else
      mean(predicted.value[is.finite(predicted.value)]), .groups='drop') |>
    pivot_wider(names_from=TRAIT, values_from=value)
  m <- as.matrix(w[,traits])
  m[!is.finite(m)] <- NA_real_
  r <- cor(m, use='pairwise.complete.obs', method='pearson')
  n <- crossprod(1L*is.finite(m))
  list(wide=w, r=r, n=n)
}
train <- prepare('blues_me_train.rds','germplasmName')
prepare_uiuc_observed <- function() {
  file <- file.path(project_dir, 'data', 'FHB_Project_Testing_Data_2025.csv')
  raw <- read.csv(file, stringsAsFactors=FALSE)
  stopifnot(all(c('ID','SUB_NURNAME','EXPT','YR',traits) %in% names(raw)))
  raw <- raw |> filter(toupper(trimws(EXPT)) == 'UIUC', YR == 2025,
                      SUB_NURNAME %in% c('Elite', 'WSS AYN'))
  if (!nrow(raw)) stop('No replicated 2025 UIUC observations found')
  if (!all(c('Elite','WSS AYN') %in% raw$SUB_NURNAME))
    stop('Expected both replicated sub-nurseries: Elite and WSS AYN')
  for (key in c('ID','SUB_NURNAME')) {
    if (any(is.na(raw[[key]]) | trimws(as.character(raw[[key]])) == ''))
      stop('Missing UIUC identifier: ', key)
  }
  if (!all(vapply(raw[traits], is.numeric, logical(1))))
    stop('UIUC trait columns must be numeric')
  w <- raw |>
    pivot_longer(all_of(traits), names_to='TRAIT', values_to='value') |>
    group_by(SUB_NURNAME, ID, TRAIT) |>
    summarise(value=if (any(is.finite(value))) mean(value[is.finite(value)]) else NA_real_,
              .groups='drop') |>
    pivot_wider(names_from=TRAIT, values_from=value)
  m <- as.matrix(w[,traits])
  m[!is.finite(m)] <- NA_real_
  list(wide=w, r=cor(m, use='pairwise.complete.obs', method='pearson'),
       n=crossprod(1L*is.finite(m)), nursery_count=n_distinct(w$SUB_NURNAME))
}
test <- prepare_uiuc_observed()
for(pop in c('train','test')) {
  x <- get(pop)
  d <- as.data.frame(as.table(x$r));names(d)<-c('Trait_1','Trait_2','Pearson_r')
  d$Paired_units <- as.vector(x$n)
  
# Pipeline-stage table folders
write.csv(d,file.path(out_dir,paste0('FHB_correlations_',pop,'_values.csv')),row.names=FALSE)
  matrix_table <- data.frame(Trait=rownames(x$r), x$r, check.names=FALSE)
  write.csv(matrix_table, file.path(out_dir,paste0('FHB_correlations_',pop,'_matrix.csv')), row.names=FALSE)
  cat('\n',pop,'\n'); print(d)
}
