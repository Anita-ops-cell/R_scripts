# script to perform standard workflow steps to analyze single cell RNA-Seq data
# data: 10k Human DTC Melanoma, Chromium GEM-X Single Cell 3'
# data source :https://www.10xgenomics.com/datasets/10k-human-dtc-melanoma-GEM-X


#set direcctory
setwd("D:\\Single cell_RNASeq\\Scripts")


# load libraries
install.packages("Seurat")
install.packages("Seurat", dependencies = TRUE)

library(Seurat)
library(tidyverse)

# Load the DTC  dataset
dtc.sparse.m <- Read10X_h5(filename = "C:\\Users\\sanga\\Downloads\\10k_Human_DTC_Melanoma_3p_gemx_10k_Human_DTC_Melanoma_3p_gemx_count_sample_filtered_feature_bc_matrix.h5")
str(dtc.sparse.m)
cts <- dtc.sparse.m

# Initialize the Seurat object with the raw (non-normalized data).
dtc.seurat.obj <- CreateSeuratObject(counts = cts, project = "HUMAN DTC MELANOMA", min.cells = 3, min.features = 200)
str(dtc.seurat.obj)

# 29104 features across 10367 samples

#  QC -------
View(dtc.seurat.obj@meta.data)


# % MT reads (percentage of mtochondrial reads found in each cell)
dtc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(dtc.seurat.obj, pattern = "^MT-")
View(dtc.seurat.obj@meta.data)


VlnPlot(dtc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(dtc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')


#  Filtering -----------------
dtc.seurat.obj <- subset(dtc.seurat.obj, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & 
                             percent.mt < 10)
str(dtc.seurat.obj)

#  Normalize data ----------
#dtc.seurat.obj <- NormalizeData(dtc.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
# OR
dtc.seurat.obj <- NormalizeData(dtc.seurat.obj)
str(dtc.seurat.obj)


# Identify higly variable feature---------
dtc.seurat.obj <- FindVariableFeatures(dtc.seurat.obj,selection.method ="vst",nfeatures = 2000)


# Identify the 10 most highly variable genes 
top10 <- head(VariableFeatures(dtc.seurat.obj), 10)
str(top10)


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(dtc.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Scaling the Data 
#Scale the data to center and normalize expression of genes.
#It prepares the data for dimensionality reduction (PCA) by scaling and centering gene expression values.
dtc.seurat.obj <- ScaleData(dtc.seurat.obj)

str(dtc.seurat.obj)


#  Perform Linear dimensionality reduction --------------
dtc.seurat.obj <- RunPCA(dtc.seurat.obj, features = VariableFeatures(object = dtc.seurat.obj))

# visualize PCA results
print(dtc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(dtc.seurat.obj, dims = 1, cells = 410, balanced = TRUE)

# determine dimensionality of the data
ElbowPlot(dtc.seurat.obj)

#  Clustering ------------
#Identify the Nearest Neighbors: Use the selected number of PCs (e.g., 15) for constructing a k-nearest neighbors (kNN) graph:
dtc.seurat.obj <- FindNeighbors(dtc.seurat.obj, dims = 1:10)

# understanding resolution
dtc.seurat.obj <- FindClusters(dtc.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1,2,3,4,5))
View(dtc.seurat.obj@meta.data)

DimPlot(dtc.seurat.obj, group.by = "RNA_snn_res.1", label = TRUE)

# setting identity of clusters
Idents(dtc.seurat.obj)
Idents(dtc.seurat.obj) <- "RNA_snn_res.1"
Idents(dtc.seurat.obj)

dtc.seurat.obj <- RunTSNE(dtc.seurat.obj, dims = 1:10)
DimPlot(dtc.seurat.obj, reduction = "tsne", label = TRUE)
# individual clusters
DimPlot(dtc.seurat.obj, reduction = "umap")

markers <- FindAllMarkers(dtc.seurat.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(markers)
# Check the number of cells per cluster
table(dtc.seurat.obj$seurat_clusters)
