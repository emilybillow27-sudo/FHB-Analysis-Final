library(data.table)
library(ggplot2)

### 1. Remove non-numeric genotype column (critical fix)
geno_numeric <- geno_mat[, -1]        # drop "genotype" column
geno_numeric <- as.data.table(geno_numeric)

### 2. Function to compute r^2 between two SNPs
calc_r2 <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) return(NA_real_)
  cor(as.numeric(x[ok]), as.numeric(y[ok]))^2
}

### 3. Compute LD (r^2) between adjacent SNPs
r2_vals <- sapply(1:(ncol(geno_numeric) - 1), function(i) {
  calc_r2(geno_numeric[[i]], geno_numeric[[i + 1]])
})

### 4. Build LD decay dataframe
ld_df <- data.frame(
  dist = seq_along(r2_vals),
  r2 = r2_vals
)

### 5. Plot LD decay
ggplot(ld_df, aes(x = dist, y = r2)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(se = FALSE, color = "blue") +
  labs(
    x = "Marker index distance",
    y = expression(r^2),
    title = "LD Decay Plot"
  ) +
  theme_minimal()