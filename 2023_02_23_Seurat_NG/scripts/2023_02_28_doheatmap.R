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

n <- 5


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
# Read in input file
##################################################

folder_path = file.path("../output/SingleR")

dat <- readRDS(file = file.path(folder_path, "data.rds"))


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

    data_markers <- tryCatch({
        read.table(
            file = file.path("../output/Markers", paste0("data_markers_", libraries[i], ".txt")),
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
            top_n_markers <- data_markers %>%
                group_by(cluster) %>%
                top_n(n = n, wt = avg_log2FC)

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
                height = 10
            )

            write.table(
                x = top_n_markers,
                file = file.path(output_path, paste0("doheatmap_", n, "_", libraries[i], ".txt")),
                sep = "\t",
                na = "",
                quote = FALSE,
                row.names = FALSE
            )
        }
    }

}
