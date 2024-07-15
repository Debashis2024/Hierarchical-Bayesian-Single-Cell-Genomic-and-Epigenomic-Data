if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Signac")
BiocManager::install("Seurat")
BiocManager::install("GenomeInfoDb")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("JASPAR2020")
BiocManager::install("TFBSTools")
install.packages("rstan")
install.packages("ggplot2")

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg19)
library(JASPAR2020)
library(TFBSTools)
library(rstan)
library(Matrix)
library(ggplot2)



set.seed(123) # for reproducibility

# Parameters
n_cells <- 500
n_peaks <- 1000

# Simulate peak counts (sparse matrix with more variability)
counts <- Matrix::rsparsematrix(n_cells, n_peaks, density = 0.2)
rownames(counts) <- paste0("Cell-", seq_len(n_cells))
colnames(counts) <- paste0("Peak-", seq_len(n_peaks))

# Ensure non-zero counts
counts@x[counts@x < 0.1] <- 0.1

# Simulate metadata
metadata <- data.frame(
  sample = rep(c("Sample1", "Sample2"), each = n_cells / 2),
  row.names = rownames(counts)
)

# Create Seurat object
pbmc <- CreateSeuratObject(
  counts = counts,
  assay = "peaks",
  meta.data = metadata
)




######################



# Normalize the data
pbmc <- NormalizeData(pbmc)

# Identify variable features
pbmc <- FindVariableFeatures(pbmc)

# Scale the data
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))

# Run PCA for dimensionality reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))

# Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP for visualization
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Plot the UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE)

# Fetch motifs from JASPAR2020
motifs <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", species = 9606)
)

# Add motif information
pbmc <- AddMotifs(
  object = pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg19,
  pfm = motifs
)

# Run chromVAR
pbmc <- RunChromVAR(
  object = pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg19
)

# Normalize the data again after adding chromVAR results
pbmc <- NormalizeData(pbmc)

# Identify variable features again
pbmc <- FindVariableFeatures(pbmc)

# Scale the data again
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))

# Run PCA again for dimensionality reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))

# Cluster the cells again
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP again for visualization
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Plot the UMAP again
DimPlot(pbmc, reduction = "umap", label = TRUE)
