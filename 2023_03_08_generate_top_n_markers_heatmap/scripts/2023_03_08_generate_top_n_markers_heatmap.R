#!/usr/bin/Rscript --vanilla
rm(list=ls())

set.seed(1)

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(argparse)

library(pheatmap)


##################################################
# Constants/Variables
##################################################


##################################################
# Argparse
##################################################

parser <- argparse::ArgumentParser()

parser$add_argument("-n", "--num", type="integer", help="Top number for slice", default=10)
parser$add_argument("-c", "--condition", type="character", help="Condition")
parser$add_argument("--remove_ensembl_ids", action='store_true', help="Remove Ensembl IDs")
parser$add_argument("--remove_mitochondria", action='store_true', help="Remove mitochondria")

args <- parser$parse_args()

condition <- tryCatch({
	args$condition
}, error = function(e) {
	return(NULL)
})

n <- tryCatch({
	args$num
}, error = function(e) {
	return(NULL)
})

remove_ensembl_ids <- tryCatch({
	args$remove_ensembl_ids
}, error = function(e) {
	return(FALSE)
})

remove_mitochondria <- tryCatch({
	args$remove_mitochondria
}, error = function(e) {
	return(FALSE)
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

folder_path = file.path(paste0("../../2023_02_23_Seurat_", condition, "/output/SignificantMarkers"))

data_markers <- read.table(
    file = file.path(folder_path, "data_markers_seurat_clusters.txt"),
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = ""
)


folder_path = file.path(paste0("../../2023_03_06_extract_normalized_counts/output"))

processed_dat <- readRDS(file = file.path(folder_path, paste0("data_", condition, ".rds")))


##################################################
# Process data
##################################################


if (remove_ensembl_ids) {
    data_markers <- data_markers %>%
        filter(!startsWith(gene, "ENS"))
}

if (remove_mitochondria) {
    data_markers <- data_markers %>%
        filter(!startsWith(gene, "MT-")) %>%
        filter(!startsWith(gene, "mt-"))
}

top_n_markers <- data_markers %>%
    group_by(cluster) %>%
    top_n(n = n, wt = avg_log2FC) %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)


mat1 <- processed_dat %>% 
    filter(Gene %in% top_n_markers$gene) %>%
    column_to_rownames(var = "Gene") %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

mat1 <- as.matrix(mat1)


##################################################
# Plotting
##################################################

paletteLength <- 50


myColor1 <- colorRampPalette(c("yellow", "red", "darkred"))(paletteLength)

myBreaks1 <- unique(
    c(
        seq(min(mat1, na.rm = TRUE), max(mat1, na.rm = TRUE), length.out=floor(paletteLength) )
    )
)


pheatmap(
    mat1,
    fontsize=6,
    color=myColor1,
    breaks=myBreaks1,
    cluster_rows=ifelse(nrow(mat1)>1, TRUE, FALSE),
    cluster_cols=ifelse(ncol(mat1)>1, TRUE, FALSE),
    fontsize_col=6,
    fontsize_row=3,
    filename=file.path(output_path, paste0("top_", n, "_markers_", condition, ".png")),
    fontsize_number = 6,
    width = 14,
    height = 7
)
