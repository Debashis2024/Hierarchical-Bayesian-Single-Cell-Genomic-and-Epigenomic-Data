# Ensure necessary packages are installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required Bioconductor packages
BiocManager::install("chromVAR")
BiocManager::install("Signac")
BiocManager::install("Seurat")
BiocManager::install("GenomeInfoDb")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("JASPAR2020")
BiocManager::install("TFBSTools")
BiocManager::install("motifmatchr")

# Install additional required packages
install.packages("rstan")
install.packages("ggplot2")
install.packages("pheatmap")

# Load necessary libraries
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
library(pheatmap)

################################
set.seed(123) # for reproducibility

# Parameters
n_cells <- 500
n_peaks <- 1000

# Simulate peak counts (sparse matrix with more variability)
counts <- Matrix::rsparsematrix(n_peaks, n_cells, density = 0.2)
rownames(counts) <- paste0("Peak-", seq_len(n_peaks))
colnames(counts) <- paste0("Cell-", seq_len(n_cells))

# Ensure non-zero counts
counts@x[counts@x < 0.1] <- 0.1

# Simulate peak ranges
peaks <- GRanges(
  seqnames = Rle(rep("chr1", n_peaks)),
  ranges = IRanges(start = seq(1, by = 2000, length.out = n_peaks), width = 100),
  strand = Rle(rep("*", n_peaks))
)

# Simulate metadata
metadata <- data.frame(
  sample = rep(c("Sample1", "Sample2"), each = n_cells / 2),
  row.names = colnames(counts)
)

# Create a ChromatinAssay object
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c("-", "-"),
  genome = 'hg19',
  ranges = peaks
)

# Create Seurat object
pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

# Normalize the data
pbmc <- NormalizeData(pbmc, assay = "peaks")

# Identify variable features
pbmc <- FindVariableFeatures(pbmc, assay = "peaks")

# Extract variable features and ensure they are aligned
variable_features <- VariableFeatures(pbmc)

# Ensure the order of features in 'data'
ordered_features <- rownames(pbmc[["peaks"]]@data)

# Align features and ensure the order matches
aligned_features <- ordered_features[ordered_features %in% variable_features]

# Reorder variable_features to match the order in aligned_features
variable_features <- variable_features[match(aligned_features, variable_features)]

# Check the alignment
if (!all(aligned_features == variable_features)) {
  stop("Features are not properly aligned.")
}

# Scale the data while aligning features
pbmc <- ScaleData(pbmc, features = variable_features, assay = "peaks")

# Run PCA for dimensionality reduction
pbmc <- RunPCA(pbmc, features = variable_features, assay = "peaks")

# Check that PCA results are added
pbmc <- ProjectDim(pbmc)

# Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP for visualization
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Plot the UMAP
umap_plot <- DimPlot(pbmc, reduction = "umap", label = TRUE)
ggsave("umap_plot.png", plot = umap_plot)

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
pbmc <- NormalizeData(pbmc, assay = "peaks")

# Identify variable features again
pbmc <- FindVariableFeatures(pbmc, assay = "peaks")

# Align features again
variable_features <- VariableFeatures(pbmc)

# Ensure the order of features in 'data'
ordered_features <- rownames(pbmc[["peaks"]]@data)

# Align features and ensure the order matches
aligned_features <- ordered_features[ordered_features %in% variable_features]

# Reorder variable_features to match the order in aligned_features
variable_features <- variable_features[match(aligned_features, variable_features)]

# Check the alignment
if (!all(aligned_features == variable_features)) {
  stop("Features are not properly aligned.")
}

# Scale the data again
pbmc <- ScaleData(pbmc, features = variable_features, assay = "peaks")

# Run PCA again for dimensionality reduction
pbmc <- RunPCA(pbmc, features = variable_features, assay = "peaks")

# Cluster the cells again
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP again for visualization
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Plot the UMAP again
umap_plot <- DimPlot(pbmc, reduction = "umap", label = TRUE)
ggsave("umap_plot_post_chromVAR.png", plot = umap_plot)

# Plot raw data (count matrix heatmap)
# Convert sparse matrix to dense for heatmap
counts_dense <- as.matrix(counts)

# Plot heatmap of raw data
heatmap_plot <- pheatmap::pheatmap(
  counts_dense,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Heatmap of Raw Count Data"
)

# Save the heatmap plot
ggsave("heatmap_plot.png", plot = heatmap_plot)
