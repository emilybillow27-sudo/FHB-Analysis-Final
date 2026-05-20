library(dplyr)
library(ggplot2)
library(stringr)

# =========================================================
# 1. STATE PROGRAM ASSIGNMENT
# =========================================================

program_lookup <- data.frame(genotype = geno_mat$genotype) %>%
  mutate(program = case_when(
    str_detect(genotype, "^CO") ~ "CO",
    str_detect(genotype, "^KS") ~ "KS",
    str_detect(genotype, "^MT") ~ "MT",
    str_detect(genotype, "^NE") ~ "NE",
    str_detect(genotype, "^OK") ~ "OK",
    str_detect(genotype, "^SD") ~ "SD",
    str_detect(genotype, "^TX") ~ "TX",
    str_detect(genotype, "^VA") ~ "VA",
    TRUE ~ "Other"
  ))

# =========================================================
# 2. TRAINING / TESTING SET ASSIGNMENT
# =========================================================

set_lookup <- bind_rows(
  test %>% distinct(ID) %>% mutate(genotype = ID, set = "Testing") %>% select(genotype, set),
  train %>% distinct(germplasmName) %>% mutate(genotype = germplasmName, set = "Training") %>% select(genotype, set)
) %>% distinct()

# =========================================================
# 3. MERGE METADATA
# =========================================================

metadata_df <- program_lookup %>%
  left_join(set_lookup, by = "genotype")

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
  left_join(metadata_df, by = "genotype")

# =========================================================
# 5. PCA PLOT WITH OKABE–ITO PALETTE
# =========================================================

ggplot(pca_df, aes(
  x = PC1,
  y = PC2,
  color = program,
  shape = set
)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_shape_manual(values = c("Training" = 16, "Testing" = 17)) +
  scale_color_manual(
    values = c(
      "CO" = "#56B4E9",   # sky blue
      "KS" = "#E69F00",   # orange
      "MT" = "#009E73",   # bluish green
      "NE" = "#D55E00",   # vermilion
      "OK" = "#CC79A7",   # reddish purple
      "Other" = "#000000",# black
      "SD" = "#F0E442",   # yellow
      "TX" = "#999999",   # gray
      "VA" = "#0072B2"    # blue
    ),
    breaks = c("CO", "KS", "MT", "NE", "OK", "Other", "SD", "TX", "VA")
  ) +
  labs(
    title = "PCA of Genotypes by Program (Color) and Set (Shape)",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )








# ====== Compute eigenvalues, variance explained, cumulative variance ======
eigs <- pca$sdev^2
var_expl <- eigs / sum(eigs) * 100
cum_var <- cumsum(var_expl)

# ====== Create PC table (first 64 PCs) ======
pc_table <- data.frame(
  PC = paste0("PC", 1:64),
  Eigenvalue = round(eigs[1:64], 3),
  PercentVariance = round(var_expl[1:64], 2),
  CumulativeVariance = round(cum_var[1:64], 2)
)

pc_table

# ====== Scree plot for PCs 1–64 (NO x-axis labels) ======
library(ggplot2)

scree_df <- data.frame(
  PC = 1:length(eigs),
  Variance = var_expl
)

# Optional: small visual confidence band (±0.1%)
scree_df$lower <- scree_df$Variance - 0.1
scree_df$upper <- scree_df$Variance + 0.1

p <- ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, fill = "steelblue") +
  geom_point(size = 2) +
  geom_line() +
  scale_x_continuous(limits = c(1, 64)) +   # show full range, no labels
  labs(
    title = "Scree Plot of Principal Components",
    x = "Principal Component",
    y = "Percent Variance Explained (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),          # remove x-axis labels
    axis.ticks.x = element_blank(),         # remove x-axis ticks
    plot.title = element_text(face = "bold")
  )

print(p)
