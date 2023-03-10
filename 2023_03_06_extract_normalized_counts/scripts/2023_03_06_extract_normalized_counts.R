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


##################################################
# Process data
##################################################

DefaultAssay(dat) <- "RNA"


processed_dat <- GetAssayData(object = dat,  slot = "data")
raw_dat <- GetAssayData(object = dat, slot = "counts")
scale_dat <- GetAssayData(object = dat, slot = "scale.data")


processed_dat <- processed_dat %>%
	as.data.frame(check.names = FALSE, stringsAsFactors = FALSE) %>%
	rownames_to_column(var = "Gene") %>%
	as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

raw_dat <- raw_dat %>%
	as.data.frame(check.names = FALSE, stringsAsFactors = FALSE) %>%
	rownames_to_column(var = "Gene") %>%
	as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

scale_dat <- scale_dat %>%
	as.data.frame(check.names = FALSE, stringsAsFactors = FALSE) %>%
	rownames_to_column(var = "Gene") %>%
	as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)


##################################################
# Save data
##################################################

saveRDS(processed_dat, file = file.path(output_path, paste0("data_", condition, ".rds")))
saveRDS(raw_dat, file = file.path(output_path, paste0("counts_", condition, ".rds")))
saveRDS(scale_dat, file = file.path(output_path, paste0("scale_data_", condition, ".rds")))
