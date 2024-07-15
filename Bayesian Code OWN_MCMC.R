# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("chromVAR", "Signac", "Seurat", "GenomeInfoDb", "BSgenome.Hsapiens.UCSC.hg19", "JASPAR2020", "TFBSTools", "motifmatchr"))
install.packages(c("rstan", "ggplot2"))

library(chromVAR)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg19)
library(JASPAR2020)
library(TFBSTools)
library(rstan)
library(Matrix)
library(ggplot2)
library(motifmatchr)

# Simulate data
set.seed(123)
n_cells <- 500
n_peaks <- 1000
counts <- Matrix::rsparsematrix(n_peaks, n_cells, density = 0.2)
rownames(counts) <- paste0("Peak-", seq_len(n_peaks))
colnames(counts) <- paste0("Cell-", seq_len(n_cells))
counts@x[counts@x < 0.1] <- 0.1
peaks <- GRanges(seqnames = Rle(rep("chr1", n_peaks)), ranges = IRanges(start = seq(1, by = 2000, length.out = n_peaks), width = 100), strand = Rle(rep("*", n_peaks)))
metadata <- data.frame(sample = rep(c("Sample1", "Sample2"), each = n_cells / 2), row.names = colnames(counts))

chrom_assay <- CreateChromatinAssay(counts = counts, sep = c("-", "-"), genome = 'hg19', ranges = peaks)
pbmc <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = metadata)

pbmc <- NormalizeData(pbmc, assay = "peaks")
pbmc <- FindVariableFeatures(pbmc, assay = "peaks")
variable_features <- VariableFeatures(pbmc)
ordered_features <- rownames(pbmc[["peaks"]]@data)
aligned_features <- ordered_features[ordered_features %in% variable_features]
variable_features <- variable_features[match(aligned_features, variable_features)]
if (!all(aligned_features == variable_features)) stop("Features are not properly aligned.")
pbmc <- ScaleData(pbmc, features = variable_features, assay = "peaks")
pbmc <- RunPCA(pbmc, features = variable_features, assay = "peaks")
pbmc <- ProjectDim(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE)
motifs <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", species = 9606))
pbmc <- AddMotifs(object = pbmc, genome = BSgenome.Hsapiens.UCSC.hg19, pfm = motifs)
pbmc <- RunChromVAR(object = pbmc, genome = BSgenome.Hsapiens.UCSC.hg19)
pbmc <- NormalizeData(pbmc, assay = "peaks")
pbmc <- FindVariableFeatures(pbmc, assay = "peaks")
variable_features <- VariableFeatures(pbmc)
ordered_features <- rownames(pbmc[["peaks"]]@data)
aligned_features <- ordered_features[ordered_features %in% variable_features]
variable_features <- variable_features[match(aligned_features, variable_features)]
if (!all(aligned_features == variable_features)) stop("Features are not properly aligned.")
pbmc <- ScaleData(pbmc, features = variable_features, assay = "peaks")
pbmc <- RunPCA(pbmc, features = variable_features, assay = "peaks")
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE)

# Custom MCMC Implementation
# Gibbs Sampler for a simplified version of the hierarchical model

# Initialize parameters
beta_0 <- rnorm(n_peaks, 0, 1)
beta_1 <- rnorm(n_peaks, 0, 1)
u <- rnorm(n_cells, 0, 1)
epsilon <- matrix(rnorm(n_peaks * n_cells, 0, 1), nrow = n_peaks, ncol = n_cells)
sigma_u <- 1
sigma_epsilon <- 1

# Gibbs Sampler function
gibbs_sampler <- function(counts, x, n_iter = 1000) {
  n_peaks <- nrow(counts)
  n_cells <- ncol(counts)
  beta_0 <- rnorm(n_peaks, 0, 1)
  beta_1 <- rnorm(n_peaks, 0, 1)
  u <- rnorm(n_cells, 0, 1)
  epsilon <- matrix(rnorm(n_peaks * n_cells, 0, 1), nrow = n_peaks, ncol = n_cells)
  sigma_u <- 1
  sigma_epsilon <- 1
  
  samples <- list(beta_0 = matrix(0, nrow = n_iter, ncol = n_peaks),
                  beta_1 = matrix(0, nrow = n_iter, ncol = n_peaks),
                  u = matrix(0, nrow = n_iter, ncol = n_cells),
                  epsilon = array(0, dim = c(n_iter, n_peaks, n_cells)),
                  sigma_u = numeric(n_iter),
                  sigma_epsilon = numeric(n_iter))
  
  for (iter in 1:n_iter) {
    # Update beta_0
    for (j in 1:n_peaks) {
      lambda_ij <- exp(beta_0[j] + beta_1[j] * x + u + epsilon[j, ])
      beta_0[j] <- rnorm(1, mean = mean(log(lambda_ij)), sd = sqrt(sigma_epsilon / n_cells))
    }
    
    # Update beta_1
    for (j in 1:n_peaks) {
      lambda_ij <- exp(beta_0[j] + beta_1[j] * x + u + epsilon[j, ])
      beta_1[j] <- rnorm(1, mean = mean(log(lambda_ij) / x), sd = sqrt(sigma_epsilon / n_cells))
    }
    
    # Update u
    for (i in 1:n_cells) {
      lambda_ij <- exp(beta_0 + beta_1 * x[i] + u[i] + epsilon[, i])
      u[i] <- rnorm(1, mean = mean(log(lambda_ij)), sd = sqrt(sigma_u / n_peaks))
    }
    
    # Update epsilon
    for (j in 1:n_peaks) {
      for (i in 1:n_cells) {
        lambda_ij <- exp(beta_0[j] + beta_1[j] * x[i] + u[i] + epsilon[j, i])
        epsilon[j, i] <- rnorm(1, mean = log(lambda_ij) - (beta_0[j] + beta_1[j] * x[i] + u[i]), sd = sigma_epsilon)
      }
    }
    
    # Update sigma_u
    sigma_u <- 1 / rgamma(1, shape = 0.5 * n_cells, rate = 0.5 * sum(u^2))
    
    # Update sigma_epsilon
    sigma_epsilon <- 1 / rgamma(1, shape = 0.5 * n_peaks * n_cells, rate = 0.5 * sum(epsilon^2))
    
    # Store samples
    samples$beta_0[iter, ] <- beta_0
    samples$beta_1[iter, ] <- beta_1
    samples$u[iter, ] <- u
    samples$epsilon[iter, , ] <- epsilon
    samples$sigma_u[iter] <- sigma_u
    samples$sigma_epsilon[iter] <- sigma_epsilon
  }
  
  return(samples)
}

# Simulate covariate data
x <- rep(1, n_cells)

# Run Gibbs sampler
samples <- gibbs_sampler(counts, x, n_iter = 1000)

# Analyze and visualize the results
beta_0_means <- colMeans(samples$beta_0)
beta_1_means <- colMeans(samples$beta_1)
u_means <-Certainly! Below is the continuation and completion of the R code with a detailed implementation of MCMC using Gibbs sampling, followed by explanations and visualizations.

### Continuation and Completion of MCMC Implementation