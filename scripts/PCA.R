library(dplyr)
library(ggplot2)
library(stringr)

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
  test %>%
    dplyr::distinct(ID) %>%
    dplyr::rename(genotype = ID) %>%
    dplyr::mutate(set = "Testing") %>%
    dplyr::select(genotype, set),
  
  train %>%
    dplyr::distinct(germplasmName) %>%
    dplyr::rename(genotype = germplasmName) %>%
    dplyr::mutate(set = "Training") %>%
    dplyr::select(genotype, set)
) %>%
  dplyr::distinct()

# =========================================================
# 3. MERGE METADATA
# =========================================================
metadata_df <- program_lookup %>%
  dplyr::left_join(set_lookup, by = "genotype")

# =========================================================
# 4. PCA
# =========================================================
geno_numeric <- geno_mat[, -1]
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
  scale_shape_manual(values = c("Training" = 16, "Testing" = 17)) +
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
