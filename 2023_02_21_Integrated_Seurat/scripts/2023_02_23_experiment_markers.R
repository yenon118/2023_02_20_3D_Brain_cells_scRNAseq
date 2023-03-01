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
library(SingleR)

library(biomaRt)


##################################################
# Constants/Variables
##################################################


##################################################
# Output folder
##################################################

output_path <- file.path("../output/ExperimentMarkers")

if(!dir.exists(output_path)){
	dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
	if(!dir.exists(output_path)){
		quit(status=1)
	}
}


##################################################
# Read in input file
##################################################

folder_path = file.path("../output/SingleR")

dat <- readRDS(file = file.path(folder_path, "data.rds"))


##################################################
# Argparse
##################################################

parser <- argparse::ArgumentParser()

parser$add_argument("-l", "--libraries", type="character", action="append", help="Libraries for annotations")

args <- parser$parse_args()


##################################################
# Process data
##################################################

DefaultAssay(dat) <- "RNA"

dat[["HumanPrimaryCellAtlasData_Experiment"]] <- paste0(dat[["HumanPrimaryCellAtlasData"]][,1], " (", dat[["Experiment"]][,1], ")")
dat[["BlueprintEncodeData_Experiment"]] <- paste0(dat[["BlueprintEncodeData"]][,1], " (", dat[["Experiment"]][,1], ")")
dat[["MouseRNAseqData_Experiment"]] <- paste0(dat[["MouseRNAseqData"]][,1], " (", dat[["Experiment"]][,1], ")")
dat[["ImmGenData_Experiment"]] <- paste0(dat[["ImmGenData"]][,1], " (", dat[["Experiment"]][,1], ")")
dat[["DatabaseImmuneCellExpressionData_Experiment"]] <- paste0(dat[["DatabaseImmuneCellExpressionData"]][,1], " (", dat[["Experiment"]][,1], ")")
dat[["NovershternHematopoieticData_Experiment"]] <- paste0(dat[["NovershternHematopoieticData"]][,1], " (", dat[["Experiment"]][,1], ")")
dat[["MonacoImmuneData_Experiment"]] <- paste0(dat[["MonacoImmuneData"]][,1], " (", dat[["Experiment"]][,1], ")")

libraries <- tryCatch({
	args$libraries
}, error = function(e) {
	return(NULL)
})

if(is.null(libraries)) {
	libraries <- c(
		"HumanPrimaryCellAtlasData_Experiment",
		"BlueprintEncodeData_Experiment",
		"MouseRNAseqData_Experiment",
		"ImmGenData_Experiment",
		"DatabaseImmuneCellExpressionData_Experiment",
		"NovershternHematopoieticData_Experiment",
		"MonacoImmuneData_Experiment"
	)
}

for (i in 1:length(libraries)) {
	Idents(dat) <- libraries[i]

	dat.markers <- FindAllMarkers(
		dat,
		min.pct = 0.25,
		logfc.threshold = 0.25
	)

	print(head(dat.markers))

	df <- dat.markers %>%
		as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

	write.table(
		x = df,
		sep = "\t",
		file = file.path(output_path, paste0("data_markers_", libraries[i], ".txt")),
		na = "",
		row.names = FALSE,
		quote = FALSE
	)
}
