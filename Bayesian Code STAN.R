library(rstan)
library(Seurat)
library(Signac)
library(ggplot2)

# Extract counts and covariates from Seurat object
N <- ncol(pbmc)  # Number of cells
G <- nrow(pbmc)  # Number of genes
y <- as.matrix(pbmc@assays$peaks@counts)

# Convert y to integers
y <- round(y)
y <- apply(y, 2, as.integer)

# Ensure y is non-negative integers
y[y < 0] <- 0

# Transpose y to have dimensions (N, G)
y <- t(y)

# Simulated true values for the baseline expression level and random effects
true_beta_0 <- rnorm(G, 0, 1)  # Simulated true values for beta_0
true_u <- rnorm(N, 0, 1)  # Simulated true values for u

# For simplicity, create a covariate matrix (all ones in this example)
x <- matrix(1, nrow = N, ncol = G)

# Prepare the data list for Stan
stan_data <- list(N = N, G = G, y = y, x = x)

# Define the Stan model
stan_model_code <- "
data {
  int<lower=0> N;            // number of cells
  int<lower=0> G;            // number of genes
  int<lower=0> y[N, G];      // observed counts
  matrix[N, G] x;            // covariates (e.g., cell type)
}
parameters {
  real beta_0[G];            // baseline expression level
  real beta_1[G];            // effect of covariates
  real<lower=0> sigma_u;     // standard deviation of random effect
  real<lower=0> sigma_e;     // standard deviation of measurement noise
  vector[N] u;               // random effect
}
model {
  beta_0 ~ normal(0, 1);
  beta_1 ~ normal(0, 1);
  u ~ normal(0, sigma_u);
  for (n in 1:N)
    for (g in 1:G)
      y[n, g] ~ poisson_log(beta_0[g] + beta_1[g] * x[n, g] + u[n]);
}
"

# Compile the model
stan_model <- stan_model(model_code = stan_model_code)

# Fit the model using MCMC sampling
fit <- sampling(stan_model, data = stan_data, iter = 1000, chains = 1)

# Extract the results
stan_results <- extract(fit)
beta_0 <- colMeans(stan_results$beta_0)
beta_1 <- colMeans(stan_results$beta_1)
sigma_u <- mean(stan_results$sigma_u)
sigma_e <- mean(stan_results$sigma_e)

# Visualize the results using ggplot2
beta_0_df <- data.frame(gene = 1:G, beta_0 = beta_0)
p1 <- ggplot(beta_0_df, aes(x = gene, y = beta_0)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Posterior Means of Beta_0", x = "Gene", y = "Beta_0")
ggsave("beta_0_posterior_means.png", plot = p1)

# Visualize random effects
u_df <- data.frame(cell = 1:N, u = rowMeans(stan_results$u))
p2 <- ggplot(u_df, aes(x = cell, y = u)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Posterior Means of Random Effects", x = "Cell", y = "Random Effect")
ggsave("random_effects_posterior_means.png", plot = p2)

# Plot the trace plots for beta_0
trace_beta_0 <- stan_results$beta_0
trace_plots <- function(trace_matrix, param_name) {
  plots <- lapply(1:ncol(trace_matrix), function(i) {
    data <- data.frame(iteration = 1:nrow(trace_matrix), value = trace_matrix[, i])
    ggplot(data, aes(x = iteration, y = value)) +
      geom_line() +
      theme_minimal() +
      labs(title = paste("Trace plot for", param_name, "[", i, "]"), x = "Iteration", y = "Value")
  })
  return(plots)
}

beta_0_trace_plots <- trace_plots(trace_beta_0, "beta_0")
for (i in 1:length(beta_0_trace_plots)) {
  ggsave(paste0("trace_plot_beta_0_", i, ".png"), plot = beta_0_trace_plots[[i]])
}

# Plot the trace plots for u
trace_u <- stan_results$u
u_trace_plots <- trace_plots(trace_u, "u")
for (i in 1:length(u_trace_plots)) {
  ggsave(paste0("trace_plot_u_", i, ".png"), plot = u_trace_plots[[i]])
}

# Compare posterior means to true values
# Plot posterior means vs. true values for beta_0
beta_0_comparison <- data.frame(gene = 1:G, posterior_mean = beta_0, true_value = true_beta_0)
p3 <- ggplot(beta_0_comparison, aes(x = true_value, y = posterior_mean)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Posterior Means vs True Values for Beta_0", x = "True Value", y = "Posterior Mean")
ggsave("beta_0_posterior_vs_true.png", plot = p3)

# Plot posterior means vs. true values for u
u_comparison <- data.frame(cell = 1:N, posterior_mean = rowMeans(stan_results$u), true_value = true_u)
p4 <- ggplot(u_comparison, aes(x = true_value, y = posterior_mean)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Posterior Means vs True Values for Random Effects", x = "True Value", y = "Posterior Mean")
ggsave("u_posterior_vs_true.png", plot = p4)

# Optional: Density plots for posterior distributions
density_plot <- function(samples, param_name) {
  data <- data.frame(value = samples)
  ggplot(data, aes(x = value)) +
    geom_density() +
    theme_minimal() +
    labs(title = paste("Density plot for", param_name), x = "Value", y = "Density")
}

# Density plot for beta_0
density_beta_0 <- lapply(1:G, function(i) density_plot(stan_results$beta_0[, i], paste("beta_0[", i, "]")))
for (i in 1:length(density_beta_0)) {
  ggsave(paste0("density_plot_beta_0_", i, ".png"), plot = density_beta_0[[i]])
}

# Density plot for u
density_u <- lapply(1:N, function(i) density_plot(stan_results$u[, i], paste("u[", i, "]")))
for (i in 1:length(density_u)) {
  ggsave(paste0("density_plot_u_", i, ".png"), plot = density_u[[i]])
}

# Save all ggplots to files
ggsave("beta_0_posterior_means.png", plot = p1)
ggsave("random_effects_posterior_means.png", plot = p2)
ggsave("beta_0_posterior_vs_true.png", plot = p3)
ggsave("u_posterior_vs_true.png", plot = p4)



###########################
# Ensure counts and metadata are aligned by checking the cell names
seurat_cell_names <- colnames(pbmc@assays$peaks@counts)
count_matrix_cell_names <- colnames(y)

# Check if cell names match and subset if necessary
if (!all(seurat_cell_names == count_matrix_cell_names)) {
  common_cells <- intersect(seurat_cell_names, count_matrix_cell_names)
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



############################
#Gene Expression Variability Analysis

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Gene expression variability
# Calculate gene expression variability
gene_variability <- apply(y, 1, var)
variability_df <- data.frame(gene = 1:G, variability = gene_variability)

# Plot gene expression variability
p_variability <- ggplot(variability_df, aes(x = gene, y = variability)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Gene Expression Variability", x = "Gene", y = "Variability")

ggsave("gene_expression_variability.png", plot = p_variability)



####################
# Verify the lengths
length(clusters)  # Length of clusters
ncol(y)  # Number of columns in y

#Cluster Specific Gene Variability

# Assuming pbmc and y are aligned initially, extract the clusters
clusters <- pbmc@meta.data$seurat_clusters
clusters <- clusters[colnames(y)]

# Check again if lengths match
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


######################################
# #Differential Expression Analysis

# Extract metadata
metadata <- pbmc@meta.data

# Differential expression between Sample1 and Sample2
sample_labels <- metadata$sample
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


########################################
#Principal Component Analysis (PCA) on Posterior Means
# Perform PCA on the posterior means of beta_0
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

ggsave("differential_expression.png", plot = p_diff_expr)

################################

#Cluster-Specific Gene Expression


# Identify genes highly expressed in specific clusters
clusters <- as.numeric(pbmc@meta.data$seurat_clusters)
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
