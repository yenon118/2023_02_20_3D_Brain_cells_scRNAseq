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


fp <- file.path(paste0("../../2023_02_23_Seurat_", condition, "/output/Markers/data_markers_seurat_clusters.txt"))

data_markers <- tryCatch({
	read.table(
		file = fp,
		sep = "\t",
		header = TRUE,
		check.names = FALSE,
		stringsAsFactors = FALSE
	)
}, error = function(e) {
	return(NULL)
})


##################################################
# Process data
##################################################

data_markers <- data_markers %>%
	filter(!startsWith(gene, "ENS")) %>%
	filter(!startsWith(gene, "MT-")) %>%
	filter(!startsWith(gene, "mt-"))


##################################################
# Plotting
##################################################

DefaultAssay(dat) <- "RNA"


top_n_markers <- data_markers %>%
	group_by(cluster) %>%
	top_n(n = 5, wt = avg_log2FC)

p <- DoHeatmap(
	dat,
	features = top_n_markers$gene,
	group.by = "seurat_clusters",
	label = TRUE
) + theme(
	legend.title = element_text(colour = "black", size=15),
	legend.text = element_text(colour = "black", size=13),
	axis.title = element_text(colour = "black"),
	axis.text = element_text(colour = "black")
)

ggsave(
	filename = paste0("doheatmap_5_seurat_clusters_", condition, ".png"),
	plot = p,
	path = output_path,
	width = 12,
	height = 8,
	units = "in",
	dpi = 300
)
