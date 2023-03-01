#!/usr/bin/Rscript --vanilla
rm(list=ls())

set.seed(1)

library(png)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)

library(Matrix)
library(Seurat)
library(patchwork)


##################################################
# Constants/Variables
##################################################

chosen_dim <- 18


##################################################
# Output folder
##################################################

output_path <- file.path("../output/Seurat")

if(!dir.exists(output_path)){
	dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
	if(!dir.exists(output_path)){
		quit(status=1)
	}
}


##################################################
# Read in input file
##################################################

folder_path = file.path("../../2023_02_20_cellranger/output/cellranger_aggr")

directory <- "NG"


##################################################
# Process data
##################################################

dat_10xobj <- Read10X(data.dir = file.path(folder_path, directory, "outs", "count", "filtered_feature_bc_matrix"))

dat <- CreateSeuratObject(counts = dat_10xobj, project = directory, min.cells = 8, min.features = 200)

dat[['Experiment']] <- directory


##################################################
# Visualize QC metrics
##################################################

dat[["percent.mt"]] <- PercentageFeatureSet(
	dat,
	pattern = "^MT-"
)

# Visualize QC metrics as a violin plot
p <- VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ggsave(
	filename = paste0("VlnPlot_init_", directory, ".png"),
	plot = p,
	path = output_path,
	width = 21,
	height = 7
)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p <- plot1 + plot2

ggsave(
	filename = paste0("FeatureScatter_init_", directory, ".png"),
	plot = p,
	path = output_path,
	width = 14,
	height = 7
)


##################################################
# Subset data
##################################################

dat <- subset(dat, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 15)


##################################################
# Visualize QC metrics after filtering
##################################################

# Visualize QC metrics as a violin plot
p <- VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ggsave(
	filename = paste0("VlnPlot_filtered_", directory, ".png"),
	plot = p,
	path = output_path,
	width = 21,
	height = 7
)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p <- plot1 + plot2

ggsave(
	filename = paste0("FeatureScatter_filtered_", directory, ".png"),
	plot = p,
	path = output_path,
	width = 14,
	height = 7
)


##################################################
# Normalize and find top features
##################################################

dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = 10000)

dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)


##################################################
# Scale data
##################################################

all.genes <- rownames(dat)

dat <- ScaleData(dat, features = all.genes)


##################################################
# Run PCA
##################################################

# Run PCA on variable features
dat <- RunPCA(dat, features = VariableFeatures(object = dat))


# Examine and visualize PCA results a few different ways
print(dat[["pca"]], dims = 1:5, nfeatures = 5)

p <- VizDimLoadings(dat, dims = 1:2, reduction = "pca")

ggsave(
	filename = "VizDimLoadings.png",
	plot = p,
	path = output_path,
	width = 14,
	height = 7
)

p <- DimPlot(dat, reduction = "pca", raster=FALSE)

ggsave(
	filename = "DimPlot_PCA.png",
	plot = p,
	path = output_path,
	width = 14,
	height = 7
)

cat(rep("\n", 2))
png(filename = file.path(output_path, "DimHeatmap_PC1.png"), width = 30, height = 10, units = "in", res = 200)
DimHeatmap(dat, dims = 1, cells = 500, balanced = TRUE)
dev.off()

cat(rep("\n", 2))
png(filename = file.path(output_path, "DimHeatmap_PC1_PC15.png"), width = 30, height = 20, units = "in", res = 200)
DimHeatmap(dat, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()


# ##################################################
# # Calculate using JackStraw and ScoreJackStraw
# ##################################################

# # NOTE: This process can take a long time for big datasets, comment out for expediency. More
# # approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# # computation time
# # dat <- JackStraw(dat, num.replicate = 500, dims = 50)
# # dat <- ScoreJackStraw(dat, dims = 1:50)

# p <- JackStrawPlot(dat, dims = 1:50)

# ggsave(
#   filename = "JackStrawPlot.png",
#   plot = p,
#   path = output_path,
#   width = 14,
#   height = 7
# )


##################################################
# Plot ElbowPlot
##################################################

p <- ElbowPlot(dat, ndims = 50)

ggsave(
	filename = "ElbowPlot.png",
	plot = p,
	path = output_path,
	width = 14,
	height = 7
)


##################################################
# Cluster the cells
##################################################

dat <- FindNeighbors(dat, dims = 1:chosen_dim)
dat <- FindClusters(dat, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
print(head(Idents(dat), 5))


##################################################
# Run non-linear dimensional reduction (UMAP/tSNE)
##################################################

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
dat <- RunUMAP(dat, dims = 1:chosen_dim, reduction = "pca")
dat <- RunTSNE(dat, dims = 1:chosen_dim, reduction = "pca")

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
p <- DimPlot(dat, reduction = "umap", label = TRUE, label.size = 6, raster=FALSE)

ggsave(
	filename = "DimPlot_UMAP.png",
	plot = p,
	path = output_path,
	width = 14,
	height = 7
)

p <- DimPlot(dat, reduction = "tsne", label = TRUE, label.size = 6, raster=FALSE)

ggsave(
	filename = "DimPlot_tSNE.png",
	plot = p,
	path = output_path,
	width = 14,
	height = 7
)

p <- DimPlot(dat, reduction = "umap", group.by = 'orig.ident', pt.size = 0.1, raster=FALSE)

ggsave(
	filename = "DimPlot_UMAP_orig_ident.png",
	plot = p,
	path = output_path,
	width = 14,
	height = 7
)

p <- DimPlot(dat, reduction = "tsne", group.by = 'orig.ident', pt.size = 0.1, raster=FALSE)

ggsave(
	filename = "DimPlot_tSNE_orig_ident.png",
	plot = p,
	path = output_path,
	width = 14,
	height = 7
)

p <- DimPlot(dat, reduction = "umap", group.by = 'Experiment', pt.size = 0.1, raster=FALSE)

ggsave(
	filename = "DimPlot_UMAP_experiment.png",
	plot = p,
	path = output_path,
	width = 14,
	height = 7
)

p <- DimPlot(dat, reduction = "tsne", group.by = 'Experiment', pt.size = 0.1, raster=FALSE)

ggsave(
	filename = "DimPlot_tSNE_experiment.png",
	plot = p,
	path = output_path,
	width = 14,
	height = 7
)

p <- DimPlot(dat, reduction = "umap", split.by = "orig.ident", raster=FALSE)

ggsave(
	filename = "DimPlot_UMAP_split_by_orig_ident.png",
	plot = p,
	path = output_path,
	width = 14,
	height = 7
)

p <- DimPlot(dat, reduction = "tsne", split.by = "orig.ident", raster=FALSE)

ggsave(
	filename = "DimPlot_tSNE_split_by_orig_ident.png",
	plot = p,
	path = output_path,
	width = 14,
	height = 7
)

p <- DimPlot(dat, reduction = "umap", split.by = "Experiment", raster=FALSE)

ggsave(
	filename = "DimPlot_UMAP_split_by_experiment.png",
	plot = p,
	path = output_path,
	width = 14,
	height = 7
)

p <- DimPlot(dat, reduction = "tsne", split.by = "Experiment", raster=FALSE)

ggsave(
	filename = "DimPlot_tSNE_split_by_experiment.png",
	plot = p,
	path = output_path,
	width = 14,
	height = 7
)


##################################################
# Save data
##################################################
saveRDS(dat, file = file.path(output_path, "data.rds"))


# ##################################################
# # Finding differentially expressed features (cluster biomarkers)
# ##################################################

# # find all markers of cluster 2
# cluster2.markers <- FindMarkers(dat, ident.1 = 2, min.pct = 0.25)
# print(head(cluster2.markers, n = 5))

# # find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers <- FindMarkers(dat, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
# print(head(cluster5.markers, n = 5))

# # find markers for every cluster compared to all remaining cells, report only the positive
# # ones
# dat.markers <- FindAllMarkers(dat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# print(head(dat.markers, n = 5))

# df <- dat.markers %>%
# 	as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

# write.table(
#   x = df,
#   sep = "\t",
#   file = file.path(output_path, "data_markers.txt"),
#   na = "",
#   row.names = FALSE,
#   quote = FALSE
# )
