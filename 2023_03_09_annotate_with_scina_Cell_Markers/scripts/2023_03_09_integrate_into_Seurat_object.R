#!/usr/bin/Rscript --vanilla
rm(list=ls())

set.seed(1)

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(argparse)

library(Matrix)
library(Seurat)
library(patchwork)


##################################################
# Constants/Variables
##################################################


##################################################
# Argparse
##################################################

parser <- argparse::ArgumentParser()

parser$add_argument("-c", "--condition", type="character", help="Condition")

args <- parser$parse_args()

condition <- tryCatch({
	args$condition
}, error = function(e) {
	return(NULL)
})


##################################################
# Output folder
##################################################

output_path <- file.path("../output")

if(!dir.exists(output_path)){
	dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
	if(!dir.exists(output_path)){
		quit(status=1)
	}
}


##################################################
# Input and reference files
##################################################

folder_path = file.path(paste0("../../2023_02_23_Seurat_", condition, "/output/SingleR"))

dat <- readRDS(file = file.path(folder_path, "data.rds"))


folder_path = file.path("../output")

ref <- read.table(
	file.path(folder_path, paste0("scina_annotation_", condition, ".txt")),
	header = TRUE,
	sep = "\t",
	check.names = FALSE,
	stringsAsFactors = FALSE
)

print(head(ref))
print(tail(ref))
print(dim(ref))


##################################################
# Process data
##################################################

DefaultAssay(dat) <- "RNA"

cells <- colnames(x = dat)

print(head(cells))

idx <- match(cells, ref[,1])
custom_annotations <- ref[idx,2]
custom_annotations[custom_annotations == ""] <- "unknown"
custom_annotations[is.null(custom_annotations)] <- "unknown"
custom_annotations[is.na(custom_annotations)] <- "unknown"


dat[["custom_annotations"]] <- custom_annotations

print(head(dat[[]]))
print(tail(dat[[]]))


plotDimPlot <- function(dat, condition, reduction = "umap", ref = "HumanPrimaryCellAtlasData") {
	p <- DimPlot(dat, reduction = reduction, label = TRUE, pt.size = 0.5, group.by = ref, raster = FALSE) + NoLegend()

	if (reduction == "umap") {
		filename = paste0("DimPlot_UMAP_singler_", condition, "_", ref, "_labeled.png")
	} else if (reduction == "tsne") {
		filename = paste0("DimPlot_tSNE_singler_", condition, "_", ref, "_labeled.png")
	}

	ggsave(
		filename = filename,
		plot = p,
		path = output_path,
		width = 12,
		height = 7,
		units = "in",
		dpi = 300
	)

	p <- DimPlot(dat, reduction = reduction, label = TRUE, pt.size = 0.5, group.by = ref, split.by = "Experiment", raster=FALSE) + NoLegend()

	if (reduction == "umap") {
		filename = paste0("DimPlot_UMAP_singler_", condition, "_", ref, "_labeled_Experiment.png")
	} else if (reduction == "tsne") {
		filename = paste0("DimPlot_tSNE_singler_", condition, "_", ref, "_labeled_Experiment.png")
	}

	ggsave(
		filename = filename,
		plot = p,
		path = output_path,
		width = 12,
		height = 7,
		units = "in",
		dpi = 300
	)

	p <- DimPlot(dat, reduction = reduction, pt.size = 0.5, group.by = ref, raster = FALSE) + 
		theme(legend.position="bottom")

	if (reduction == "umap") {
		filename = paste0("DimPlot_UMAP_singler_", condition, "_", ref, ".png")
	} else if (reduction == "tsne") {
		filename = paste0("DimPlot_tSNE_singler_", condition, "_", ref, ".png")
	}

	ggsave(
		filename = filename,
		plot = p,
		path = output_path,
		width = 12,
		height = 7,
		units = "in",
		dpi = 300
	)

	p <- DimPlot(dat, reduction = reduction, pt.size = 0.5, group.by = ref, split.by = "Experiment", raster=FALSE) + 
		theme(legend.position="bottom")

	if (reduction == "umap") {
		filename = paste0("DimPlot_UMAP_singler_", condition, "_", ref, "_Experiment.png")
	} else if (reduction == "tsne") {
		filename = paste0("DimPlot_tSNE_singler_", condition, "_", ref, "_Experiment.png")
	}

	ggsave(
		filename = filename,
		plot = p,
		path = output_path,
		width = 12,
		height = 7,
		units = "in",
		dpi = 300
	)
}


plotDimPlot(dat, condition, "umap", "custom_annotations")
plotDimPlot(dat, condition, "tsne", "custom_annotations")


##################################################
# Save data
##################################################
saveRDS(dat, file = file.path(output_path, paste0("data_", condition, ".rds")))
