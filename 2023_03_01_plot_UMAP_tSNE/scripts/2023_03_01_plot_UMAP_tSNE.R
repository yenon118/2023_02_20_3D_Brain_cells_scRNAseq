#!/usr/bin/Rscript --vanilla
rm(list=ls())

set.seed(1)

library(png)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(argparse)

library(Matrix)
library(Seurat)
library(patchwork)
library(celldex)
library(SingleR)

library(openxlsx)


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
# Read in input file
##################################################

dat <- readRDS(
    file = file.path(paste0("../../2023_02_23_Seurat_", condition, "/output/SingleR/data.rds"))
)


##################################################
# Plotting
##################################################

DefaultAssay(dat) <- "RNA"


p <- DimPlot(dat, reduction = "umap", label = TRUE, label.size = 8, raster=FALSE) + 
theme(
	legend.title = element_text(colour = "black", size=15),
	legend.text = element_text(colour = "black", size=15),
	axis.title = element_text(colour = "black", size=15),
	axis.text = element_text(colour = "black", size=15)
)

ggsave(
	filename = paste0("DimPlot_UMAP_", condition, ".png"),
	plot = p,
	path = output_path,
	width = 10,
	height = 5,
    units = "in",
    dpi = 300
)

p <- DimPlot(dat, reduction = "tsne", label = TRUE, label.size = 8, raster=FALSE) + 
theme(
	legend.title = element_text(colour = "black", size=15),
	legend.text = element_text(colour = "black", size=15),
	axis.title = element_text(colour = "black", size=15),
	axis.text = element_text(colour = "black", size=15)
)

ggsave(
	filename = paste0("DimPlot_tSNE_", condition, ".png"),
	plot = p,
	path = output_path,
	width = 10,
	height = 5,
    units = "in",
    dpi = 300
)
