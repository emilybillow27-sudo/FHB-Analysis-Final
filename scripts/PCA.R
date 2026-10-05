# Run from the project root (FHB Analysis Final).
if (!file.exists("data/FHB_Project_Training_Data.csv")) {
  stop("Set the working directory to the FHB Analysis Final project root.")
}
for (output_dir in c("results/intermediate", "results/tables", "results/figures")) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

library(dplyr)
library(ggplot2)
library(stringr)

# The main pipeline saves the marker matrix used for its GRM.
if (!file.exists("results/intermediate/geno_mat.rds")) {
  stop("Missing geno_mat.rds. Run the marker preparation section of working_script_fixed.R first.")
}
geno_mat <- readRDS("results/intermediate/geno_mat.rds")
train <- readRDS("results/intermediate/blues_me_train.rds")
test <- readRDS("results/intermediate/blues_se_test.rds")

# =========================================================
# 0. FIX GENOTYPE COLUMN NAME IN geno_mat
# =========================================================
# Use the first column of geno_mat as the genotype ID (you showed it's FullSampleName)
geno_id_col <- colnames(geno_mat)[1]
geno_mat <- geno_mat %>% dplyr::rename(genotype = dplyr::all_of(geno_id_col))

# =========================================================
# 1. STATE PROGRAM ASSIGNMENT
# =========================================================
program_lookup <- geno_mat %>%
  dplyr::distinct(genotype) %>%
  dplyr::mutate(program = dplyr::case_when(
    stringr::str_detect(genotype, "^CO") ~ "CO",
    stringr::str_detect(genotype, "^KS") ~ "KS",
    stringr::str_detect(genotype, "^MT") ~ "MT",
    stringr::str_detect(genotype, "^NE") ~ "NE",
    stringr::str_detect(genotype, "^OK") ~ "OK",
    stringr::str_detect(genotype, "^SD") ~ "SD",
    stringr::str_detect(genotype, "^TX") ~ "TX",
    stringr::str_detect(genotype, "^VA") ~ "VA",
    TRUE ~ "Other"
  ))

# =========================================================
# 2. TRAINING / TESTING SET ASSIGNMENT
# =========================================================
set_lookup <- dplyr::bind_rows(
  test %>% dplyr::filter(!is.na(FullSampleName)) %>%
    dplyr::distinct(FullSampleName) %>%
    dplyr::transmute(genotype = FullSampleName, set = "Testing"),
  train %>% dplyr::filter(!is.na(FullSampleName)) %>%
    dplyr::distinct(FullSampleName) %>%
    dplyr::transmute(genotype = FullSampleName, set = "Training")
) %>%
  dplyr::group_by(genotype) %>%
  dplyr::summarise(set = paste(sort(unique(set)), collapse = " and "), .groups = "drop")

# =========================================================
# 3. MERGE METADATA
# =========================================================
metadata_df <- program_lookup %>%
  dplyr::left_join(set_lookup, by = "genotype") %>%
  dplyr::mutate(set = dplyr::coalesce(set, "Other genotyped"))

# =========================================================
# 4. PCA
# =========================================================
geno_numeric <- as.matrix(geno_mat[, -1, drop = FALSE])
if (anyNA(geno_numeric)) stop("PCA marker input contains missing values.")
geno_numeric <- geno_numeric[, apply(geno_numeric, 2, sd) > 0, drop = FALSE]
pca <- prcomp(geno_numeric, scale. = TRUE)

var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
pc1_var <- round(var_expl[1] * 100, 1)
pc2_var <- round(var_expl[2] * 100, 1)

pca_df <- data.frame(
  genotype = geno_mat$genotype,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
) %>%
  dplyr::left_join(metadata_df, by = "genotype")

# =========================================================
# 5. PCA PLOT
# =========================================================
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = program, shape = set)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_shape_manual(values = c("Training" = 16, "Testing" = 17, "Testing and Training" = 15, "Other genotyped" = 3)) +
  scale_color_manual(values = c(
    "CO"="#56B4E9","KS"="#E69F00","MT"="#009E73","NE"="#D55E00",
    "OK"="#CC79A7","Other"="#000000","SD"="#F0E442","TX"="#999999","VA"="#0072B2"
  )) +
  labs(
    title = "PCA of Genotypes by Program and Set",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  ) +
  theme_minimal(base_size = 14)

ggsave("results/figures/PCA_plot.png", p_pca, width = 8, height = 6, dpi = 300)

# =========================================================
# 6. SCREE PLOT
# =========================================================
eigs <- pca$sdev^2
var_expl <- eigs / sum(eigs) * 100
cum_var <- cumsum(var_expl)

scree_df <- data.frame(
  PC = 1:length(eigs),
  Variance = var_expl,
  lower = pmax(var_expl - 0.1, 0),
  upper = var_expl + 0.1
) %>%
  dplyr::filter(PC <= 64)   # <‑‑ THIS FIXES THE WARNINGS

p_scree <- ggplot(scree_df, aes(PC, Variance)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.15) +
  geom_point(size = 2) +
  geom_line() +
  labs(
    title = "Scree Plot of Principal Components",
    x = "Principal Component",
    y = "Percent Variance Explained (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave("results/figures/PCA_scree_plot.png", p_scree, width = 8, height = 6, dpi = 300)

write.csv(pca_df, "results/tables/PCA_scores.csv", row.names = FALSE)
write.csv(scree_df, "results/tables/PCA_variance_explained.csv", row.names = FALSE)
