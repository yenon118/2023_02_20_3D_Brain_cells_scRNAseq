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

output_path <- file.path(paste0("../output/SeparatedDimPlot/", condition))

if(!dir.exists(output_path)){
	dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
	if(!dir.exists(output_path)){
		quit(status=1)
	}
}


##################################################
# Read in input file
##################################################

folder_path <- file.path("../../2023_03_09_annotate_with_scina_Cell_Markers/output")

dat <- readRDS(file = file.path(folder_path, paste0("data_", condition, ".rds")))


##################################################
# Process data
##################################################

DefaultAssay(dat) <- "RNA"

major_group <- dat[["custom_annotations"]][,1]

unique_major_group <- unique(sort(major_group))

print(unique_major_group)


for (i in 1:length(unique_major_group)) {
	unique_major_group_dat <- subset(dat, subset = custom_annotations == unique_major_group[i])


	if (nrow(unique_major_group_dat[[]]) > 1) {
		p <- DimPlot(unique_major_group_dat, reduction = "umap", pt.size = 0.5, group.by = "custom_annotations", split.by = "Experiment", raster=FALSE) +
			theme(legend.position="bottom")

		ggsave(
			filename = paste0("DimPlot_UMAP_", condition, "_", gsub(" ", "_", unique_major_group[i]), ".png"),
			plot = p,
			path = output_path,
			width = 12,
			height = 7,
			units = "in",
			dpi = 300
		)

		p <- DimPlot(unique_major_group_dat, reduction = "tsne", pt.size = 0.5, group.by = "custom_annotations", split.by = "Experiment", raster=FALSE) +
			theme(legend.position="bottom")

		ggsave(
			filename = paste0("DimPlot_tSNE_", condition, "_", gsub(" ", "_", unique_major_group[i]), ".png"),
			plot = p,
			path = output_path,
			width = 12,
			height = 7,
			units = "in",
			dpi = 300
		)
	}
  	
}

