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

library(biomaRt)


##################################################
# Constants/Variables
##################################################


##################################################
# Argparse
##################################################

parser <- argparse::ArgumentParser()

parser$add_argument("-n", "--num", type="integer", help="Top number for slice", default=10)
parser$add_argument("-o", "--organ_type", type="character", help="Organ type", default="all")
parser$add_argument("--remove_ensembl_ids", action='store_true', help="Remove Ensembl IDs")

args <- parser$parse_args()

n <- tryCatch({
	args$num
}, error = function(e) {
	return(NULL)
})

organ_type <- tryCatch({
	args$organ_type
}, error = function(e) {
	return(NULL)
})

remove_ensembl_ids <- tryCatch({
	args$remove_ensembl_ids
}, error = function(e) {
	return(FALSE)
})


##################################################
# Output folder
##################################################

output_path <- file.path("../output/DoHeatmap")

if(!dir.exists(output_path)){
	dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
	if(!dir.exists(output_path)){
		quit(status=1)
	}
}


##################################################
# Query annotations from biomaRt
##################################################

getBMTable <- function(values, dataset) {
	ensembl <- useEnsembl(biomart = "genes", dataset = dataset)

	attributes = listAttributes(ensembl)

	# head(searchAttributes(mart = ensembl, pattern = "id"))
	# head(searchFilters(mart = ensembl, pattern = "id"))

	bm_table <- getBM(
		attributes = c(
			'ensembl_gene_id',
			'external_gene_name',
			'chromosome_name',
			'start_position',
			'end_position',
			'strand',
			'gene_biotype',
			'description'
		),
		filters = 'external_gene_name',
		values = values,
		mart = ensembl
	)

	return(bm_table)
}


##################################################
# Read in input file
##################################################

folder_path = file.path("../output/SingleR")

dat <- readRDS(file = file.path(folder_path, "data.rds"))

ref <- tryCatch({
	read.table(
		file = file.path("../data", "PanglaoDB_markers_27_Mar_2020.tsv"),
		sep = "\t",
		header = TRUE,
		check.names = FALSE,
		stringsAsFactors = FALSE
	)
}, error = function(e) {
	return(NULL)
})

if (!is.null(ref)) {
	if (organ_type != "all") {
		ref <- ref[ref$organ == organ_type, ]
	}
}


##################################################
# Process data
##################################################

DefaultAssay(dat) <- "RNA"


libraries <- c(
	"seurat_clusters",
	"HumanPrimaryCellAtlasData",
	"BlueprintEncodeData",
	"MouseRNAseqData",
	"ImmGenData",
	"DatabaseImmuneCellExpressionData",
	"NovershternHematopoieticData",
	"MonacoImmuneData",
	"Experiment"
)


for (i in 1:length(libraries)) {

	fp <- file.path("../output/Markers", paste0("data_markers_", libraries[i], ".txt"))

	if (file.exists(fp)) {

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

		if(!is.null(data_markers)) {
			if(nrow(data_markers) > 0 & ncol(data_markers) > 0) {

				if (remove_ensembl_ids) {
					top_n_markers <- data_markers %>%
						filter(!startsWith(gene, "ENS")) %>%
						group_by(cluster) %>%
						top_n(n = n, wt = avg_log2FC)
				} else {
					top_n_markers <- data_markers %>%
						group_by(cluster) %>%
						top_n(n = n, wt = avg_log2FC)
				}

				p <- DoHeatmap(
					dat,
					features = top_n_markers$gene,
					group.by = libraries[i],
					label = ifelse(libraries[i] == "seurat_clusters", TRUE, FALSE)
				)

				ggsave(
					filename = paste0("doheatmap_", n, "_", libraries[i], ".png"),
					plot = p,
					path = output_path,
					width = 14,
					height = ifelse(n > 10, 18, 10)
				)

				write.table(
					x = top_n_markers,
					file = file.path(output_path, paste0("doheatmap_", n, "_", libraries[i], ".txt")),
					sep = "\t",
					na = "",
					quote = FALSE,
					row.names = FALSE
				)

				bm_table <- tryCatch({
					getBMTable(values = unique(sort(top_n_markers$gene)), dataset = "hsapiens_gene_ensembl")
				}, error=function(e) {
					message(e)
					return(NULL)
				})

				if (!is.null(bm_table)) {
					top_n_markers$external_gene_name <- top_n_markers$gene

					annotated_top_n_markers <- top_n_markers %>%
						left_join(bm_table, by = "external_gene_name") %>%
						as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

					if (!is.null(ref)) {
						if(nrow(ref) > 0 & ncol(ref) > 0) {
							annotated_top_n_markers <- annotated_top_n_markers %>%
								left_join(ref, by = c("external_gene_name" = "official gene symbol")) %>%
								as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)
						}
					}

					write.table(
						x = annotated_top_n_markers,
						file = file.path(output_path, paste0("doheatmap_", n, "_", libraries[i], "_annotated.txt")),
						sep = "\t",
						na = "",
						quote = FALSE,
						row.names = FALSE
					)
				} else {

					if (!is.null(ref)) {
						if(nrow(ref) > 0 & ncol(ref) > 0) {
							annotated_top_n_markers <- annotated_top_n_markers %>%
								left_join(ref, by = c("external_gene_name" = "official gene symbol")) %>%
								as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

							write.table(
								x = annotated_top_n_markers,
								file = file.path(output_path, paste0("doheatmap_", n, "_", libraries[i], "_annotated.txt")),
								sep = "\t",
								na = "",
								quote = FALSE,
								row.names = FALSE
							)
						}
					}

				}

			}
		}

	}

}
