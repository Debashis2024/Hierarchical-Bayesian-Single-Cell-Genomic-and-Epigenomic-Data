# Load necessary libraries
library(Signac)
library(Seurat)
library(BiocManager)
BiocManager::install("GenomeInfoDb")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("JASPAR2020")
BiocManager::install("TFBSTools")

library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg19)
library(JASPAR2020)
library(TFBSTools)
library(rstan)
library(Matrix)
library(ggplot2)

# Load example data from Signac package
data("pbmc_scATAC")

# Create a Seurat object
seurat_obj <- CreateSeuratObject(
  counts = pbmc_scATAC@assays$peaks@counts,
  assay = "peaks"
)

# Add chromatin accessibility assay
seurat_obj <- AddChromatinAssay(
  object = seurat_obj,
  counts = pbmc_scATAC@assays$peaks@counts,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = './fragments.tsv.gz'
)

# Fetch motifs from JASPAR2020
motifs <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", species = 9606)
)

# Add motif information
seurat_obj <- AddMotifs(
  object = seurat_obj,
  genome = BSgenome.Hsapiens.UCSC.hg19,
  pfm = motifs
)

# Run chromVAR
seurat_obj <- RunChromVAR(
  object = seurat_obj,
  genome = BSgenome.Hsapiens.UCSC.hg19
)

# Normalize the data
seurat_obj <- NormalizeData(seurat_obj)

# Identify variable features
seurat_obj <- FindVariableFeatures(seurat_obj)

# Scale the data
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

# Run PCA for dimensionality reduction
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

# Cluster the cells
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Run UMAP for visualization
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Plot the UMAP
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

# Bayesian hierarchical model using rstan
# Define the model in Stan language
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

# Compile the Stan model
stan_model <- stan_model(model_code = stan_model_code)

# Prepare data for Stan
N <- ncol(seurat_obj)
G <- nrow(seurat_obj)
y <- as.matrix(seurat_obj@assays$peaks@counts)
x <- as.matrix(seurat_obj@meta.data) # assuming covariates are in meta data

# Run MCMC sampling
stan_data <- list(N = N, G = G, y = y, x = x)
fit <- sampling(stan_model, data = stan_data, iter = 1000, chains = 4)

# Summarize the results
print(fit)

# Plot the results
stan_results <- extract(fit)
beta_0 <- colMeans(stan_results$beta_0)
beta_1 <- colMeans(stan_results$beta_1)
sigma_u <- mean(stan_results$sigma_u)
sigma_e <- mean(stan_results$sigma_e)

# Visualize results using ggplot2
beta_0_df <- data.frame(gene = 1:G, beta_0 = beta_0)
ggplot(beta_0_df, aes(x = gene, y = beta_0)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Posterior Means of Beta_0", x = "Gene", y = "Beta_0")
