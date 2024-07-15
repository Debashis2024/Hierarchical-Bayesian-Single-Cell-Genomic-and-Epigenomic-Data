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


# Download dataset from GEO
# Install and load GEOquery package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")
library(GEOquery)

# Download the supplementary files for GSE129785
getGEOSuppFiles("GSE129785")

# Unzip the downloaded files
untar("GSE129785/GSE129785_RAW.tar", exdir = "GSE129785")

getGEOSuppFiles("GSE129785")

# Unzip the downloaded files
untar("GSE129785/GSE129785_RAW.tar", exdir = "GSE129785")

# Load the data
counts <- Read10X_h5("GSE129785/GSM3722718_atac_pbmc_5k_filtered_peak_bc_matrix.h5")
metadata <- read.csv("GSE129785/GSM3722718_atac_pbmc_5k_singlecell.csv", header = TRUE, row.names = 1)

# Create Seurat object
pbmc <- CreateSeuratObject(
  counts = counts,
  assay = "peaks",
  meta.data = metadata
)

# Add fragment information
pbmc <- SetFragments(
  object = pbmc,
  file = "GSE129785/GSM3722718_atac_pbmc_5k_fragments.tsv.gz"
)
