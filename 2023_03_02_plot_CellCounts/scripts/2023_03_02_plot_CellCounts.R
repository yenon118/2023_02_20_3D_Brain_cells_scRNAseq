#!/usr/bin/Rscript --vanilla
rm(list=ls())

set.seed(1)

library(png)
library(scales)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(patchwork)

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

folder_path <- file.path(paste0("../../2023_02_23_Seurat_", condition, "/output/CellCounts"))

dat_file_path <- file.path(folder_path, "cell_counts.xlsx")

dat_sheet_names <- getSheetNames(dat_file_path)

print(dat_sheet_names)


##################################################
# Process data
##################################################

plotStackBarPlot <- function(file_path_1, sheet, condition, is.num = FALSE) {

    dat1 <- read.xlsx(xlsxFile = file_path_1, sheet = sheet, skipEmptyRows = TRUE)

    colnames(dat1)[1] <- "Experiment"
    colnames(dat1)[2] <- "Cluster"

    dat1 <- dat1 %>%
        group_by(Experiment) %>%
        mutate(Percentage = round(100*Cell_Count/sum(Cell_Count, na.rm = TRUE), digits=2)) %>%
        as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

    if (is.num) {
        dat1$Cluster <- factor(dat1$Cluster, levels = sort(unique(as.numeric(as.character(dat1$Cluster)))))
    }

    if (max(dat1$Cell_Count, na.rm = TRUE) < 10) {
        by = 1
    } else if (max(dat1$Cell_Count, na.rm = TRUE) < 100) {
        by = 10
    } else if (max(dat1$Cell_Count, na.rm = TRUE) < 1000) {
        by = 50
    } else if (max(dat1$Cell_Count, na.rm = TRUE) < 5000) {
        by = 100
    } else if (max(dat1$Cell_Count, na.rm = TRUE) < 10000) {
        by = 200
    } else if (max(dat1$Cell_Count, na.rm = TRUE) < 25000) {
        by = 500
    } else if (max(dat1$Cell_Count, na.rm = TRUE) < 50000) {
        by = 1000
    } else if (max(dat1$Cell_Count, na.rm = TRUE) < 100000) {
        by = 1500
    }

    p1 <- ggplot(dat1, aes(y=Cell_Count, x=Cluster)) +
        geom_bar(position=position_dodge(width = 0.9), stat="identity") +
        geom_text(aes(label = Cell_Count), position=position_dodge(width = 0.9), vjust=-0.6, stat="identity", size=3) +
        labs(title=paste0(condition)) +
        scale_y_continuous(breaks = seq(0, max(dat1$Cell_Count, na.rm = TRUE), by=by)) +
        theme(axis.text.x = element_text(angle=45, hjust=1))

    p2 <- ggplot(dat1, aes(y=Percentage, x=Cluster)) +
        geom_bar(position=position_dodge(width = 0.9), stat="identity") +
        geom_text(aes(label=Percentage), position=position_dodge(width = 0.9), vjust=-0.6, stat="identity", size=3) +
        scale_y_continuous(breaks = seq(0, max(dat1$Percentage, na.rm = TRUE), by=10)) +
        theme(axis.text.x = element_text(angle=45, hjust=1))

    p <- p1 + p2

    ggsave(
        filename = paste0("cell_counts_", sheet, "_", condition, ".png"),
        plot = p,
        path = output_path,
        width = 10,
        height = 7,
        units = "in",
        dpi = 300
    )
}

plotStackBarPlot(dat_file_path, "seurat_clusters", condition, is.num = TRUE)
