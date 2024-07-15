# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("chromVAR", "Signac", "Seurat", "GenomeInfoDb", "BSgenome.Hsapiens.UCSC.hg19", "JASPAR2020", "TFBSTools", "motifmatchr"))
install.packages(c("rstan", "ggplot2"))

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

# Download and extract the dataset (example URL, please replace with actual URL)
url <- "https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_3k_nextgem/atac_pbmc_3k_nextgem_filtered_peak_bc_matrix.h5"
destfile <- "atac_pbmc_3k_nextgem_filtered_peak_bc_matrix.h5"
download.file(url, destfile, mode="wb")

# Load the dataset into Seurat
counts <- Read10X_h5(destfile)
metadata <- data.frame(row.names = colnames(counts), sample = "PBMC")

# Create ChromatinAssay
chrom_assay <- CreateChromatinAssay(counts = counts, sep = c("-", "-"), genome = 'hg19')
pbmc <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = metadata)

# Normalize the data
pbmc <- NormalizeData(pbmc, assay = "peaks")

# Identify variable features
pbmc <- FindVariableFeatures(pbmc, assay = "peaks")

# Scale the data
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc), assay = "peaks")

# Run PCA for dimensionality reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), assay = "peaks")

# Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP for visualization
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Plot the UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE)

# Fetch motifs from JASPAR2020
motifs <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", species = 9606))

# Add motif information
pbmc <- AddMotifs(object = pbmc, genome = BSgenome.Hsapiens.UCSC.hg19, pfm = motifs)

# Run chromVAR
pbmc <- RunChromVAR(object = pbmc, genome = BSgenome.Hsapiens.UCSC.hg19)

# Normalize the data again after adding chromVAR results
pbmc <- NormalizeData(pbmc, assay = "peaks")

# Identify variable features again
pbmc <- FindVariableFeatures(pbmc, assay = "peaks")

# Scale the data again
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc), assay = "peaks")

# Run PCA again for dimensionality reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), assay = "peaks")

# Cluster the cells again
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP again for visualization
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Plot the UMAP again
DimPlot(pbmc, reduction = "umap", label = TRUE)
