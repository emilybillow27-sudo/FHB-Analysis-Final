library(tidyverse)
library(reshape2)

# --- 1. Cluster the GRM ---
hc <- hclust(as.dist(1 - GRM))
GRM_ordered <- GRM[hc$order, hc$order]

# --- 2. Convert to long format ---
grm_long <- melt(GRM_ordered)
colnames(grm_long) <- c("Line1", "Line2", "Relationship")

# --- 3. Plot clustered heatmap ---
p <- ggplot(grm_long, aes(x = Line1, y = Line2, fill = Relationship)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, limits = c(-1, 1),
    name = "Genomic Relationship"
  ) +
  labs(
    title = "Clustered Genomic Relationship Matrix Heatmap",
    x = "Lines (clustered)",
    y = "Lines (clustered)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

# --- 4. Save figure ---
ggsave(
  filename = "results/figures/GRM_heatmap_clustered.png",
  plot = p,
  width = 8,
  height = 7,
  dpi = 300
)
