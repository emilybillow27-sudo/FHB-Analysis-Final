library(data.table)
library(ggplot2)
library(VariantAnnotation)
library(dplyr)

### 0. Build marker map from VCF and align to geno_mat

# Map from VCF
map_vcf <- data.frame(
  marker = rownames(vcf),
  CHR    = as.character(seqnames(rowRanges(vcf))),
  POS    = as.numeric(start(rowRanges(vcf))),
  stringsAsFactors = FALSE
)

# Markers used in genotype matrix (drop genotype column)
geno_markers <- colnames(geno_mat)[-1]

# Keep only markers present in both
map_use <- map_vcf[map_vcf$marker %in% geno_markers, ]

# Reorder map to match geno_mat column order
map_use <- map_use[match(geno_markers, map_use$marker), ]

# Sanity check: same length and no NA positions
stopifnot(length(geno_markers) == nrow(map_use))
stopifnot(!any(is.na(map_use$POS)))

### 1. Numeric genotype matrix (samples x markers)

geno_numeric <- as.data.table(geno_mat[, -1])

### 2. r^2 function

calc_r2 <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) return(NA_real_)
  cor(as.numeric(x[ok]), as.numeric(y[ok]))^2
}

### 3. Compute LD between adjacent SNPs on the same chromosome

# Distance between adjacent markers; NA when chromosome changes
map_use <- map_use |>
  arrange(CHR, POS) |>
  group_by(CHR) |>
  mutate(
    idx       = row_number(),                    # index within chromosome
    dist_bp   = POS - lag(POS),                  # distance to previous marker
    prev_mark = lag(marker)
  ) |>
  ungroup()

# Keep rows where there is a previous marker on same chromosome
map_pairs <- map_use |> filter(!is.na(dist_bp) & dist_bp >= 0)

# Match these marker pairs to geno_numeric columns
idx_curr <- match(map_pairs$marker,     geno_markers)
idx_prev <- match(map_pairs$prev_mark,  geno_markers)

r2_vals <- mapply(function(i, j) {
  calc_r2(geno_numeric[[j]], geno_numeric[[i]])
}, idx_curr, idx_prev)

ld_df <- data.frame(
  dist_bp = map_pairs$dist_bp,
  dist_kb = map_pairs$dist_bp / 1000,
  r2      = r2_vals
)

### 4. LD decay plot (physical distance)

p_ld <- ggplot(ld_df, aes(x = dist_kb, y = r2)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(se = FALSE, color = "blue") +
  labs(
    x = "Physical distance between adjacent SNPs (kb)",
    y = expression(r^2),
    title = "LD decay as a function of physical distance"
  ) +
  theme_minimal()

p_ld

### 5. (Optional) Save table and figure

if (!dir.exists("results/figures")) dir.create("results/figures", recursive = TRUE)
if (!dir.exists("results/tables"))  dir.create("results/tables",  recursive = TRUE)

write.csv(ld_df, "results/tables/ld_decay_physical_adjacent.csv", row.names = FALSE)

ggsave("results/figures/ld_decay_physical_adjacent.png", p_ld,
       width = 8, height = 5, dpi = 300)