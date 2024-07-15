# Load necessary libraries
library(ggplot2)
library(dplyr)

# Ensure counts and metadata are aligned by checking the cell names
seurat_cell_names <- colnames(pbmc@assays$peaks@counts)
count_matrix_cell_names <- colnames(y)

# Check if cell names match and subset if necessary
common_cells <- intersect(seurat_cell_names, count_matrix_cell_names)
if (length(common_cells) < length(seurat_cell_names) || length(common_cells) < length(count_matrix_cell_names)) {
  pbmc <- subset(pbmc, cells = common_cells)
  y <- y[, common_cells]
  clusters <- pbmc@meta.data$seurat_clusters
}

# Verify the lengths again
if (length(clusters) != ncol(y)) {
  stop("Length of clusters and number of columns in y do not match.")
}

# Identify genes highly expressed in specific clusters
cluster_specific_genes <- apply(y, 1, function(counts) {
  aov_res <- aov(counts ~ as.factor(clusters))
  summary(aov_res)[[1]][["Pr(>F)"]][1]
})

# Adjust for multiple testing (Bonferroni correction)
adjusted_cluster_p_values <- p.adjust(cluster_specific_genes, method = "bonferroni")
cluster_specific_genes <- which(adjusted_cluster_p_values < 0.05)

# Create a dataframe for visualization
cluster_specific_df <- data.frame(
  gene = cluster_specific_genes,
  p_value = adjusted_cluster_p_values[cluster_specific_genes]
)

# Plot cluster-specific genes
p_cluster_specific <- ggplot(cluster_specific_df, aes(x = gene, y = -log10(p_value))) +
  geom_point() +
  theme_minimal() +
  labs(title = "Cluster-Specific Gene Expression", x = "Gene", y = "-log10(Adjusted p-value)")

ggsave("cluster_specific_expression.png", plot = p_cluster_specific)

# Gene expression variability
gene_variability <- apply(y, 1, var)
variability_df <- data.frame(gene = 1:G, variability = gene_variability, row.names = NULL)

# Plot gene expression variability
p_variability <- ggplot(variability_df, aes(x = gene, y = variability)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Gene Expression Variability", x = "Gene", y = "Variability")

ggsave("gene_expression_variability.png", plot = p_variability)

# Differential expression analysis
sample_labels <- pbmc@meta.data$sample
differential_expression <- apply(y, 1, function(counts) {
  t.test(counts[sample_labels == "Sample1"], counts[sample_labels == "Sample2"])$p.value
})

# Adjust for multiple testing (Bonferroni correction)
adjusted_p_values <- p.adjust(differential_expression, method = "bonferroni")
diff_expr_genes <- which(adjusted_p_values < 0.05)

# Create a dataframe for visualization
diff_expr_df <- data.frame(
  gene = diff_expr_genes,
  p_value = adjusted_p_values[diff_expr_genes]
)

# Plot differentially expressed genes
p_diff_expr <- ggplot(diff_expr_df, aes(x = gene, y = -log10(p_value))) +
  geom_point() +
  theme_minimal() +
  labs(title = "Differential Expression Analysis", x = "Gene", y = "-log10(Adjusted p-value)")

ggsave("differential_expression.png", plot = p_diff_expr)

# Principal Component Analysis (PCA) on Posterior Means
pca <- prcomp(t(stan_results$beta_0), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  sample = sample_labels
)

# Plot PCA results
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = sample)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA on Posterior Means of Beta_0", x = "PC1", y = "PC2")

ggsave("pca_posterior_means.png", plot = p_pca)

