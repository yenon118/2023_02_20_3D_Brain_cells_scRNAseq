#!/usr/bin/Rscript --vanilla
rm(list=ls())

set.seed(1)

library(png)
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
# Output folder
##################################################

output_path <- file.path("../output/CellCountsStackedBars")

if(!dir.exists(output_path)){
	dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
	if(!dir.exists(output_path)){
		quit(status=1)
	}
}


##################################################
# Read in input file
##################################################

folder_path <- file.path("../output/CellCounts")

dat_file_path <- file.path(folder_path, "cell_counts.xlsx")

dat_sheet_names <- getSheetNames(dat_file_path)

print(dat_sheet_names)


##################################################
# Process data
##################################################

plotStackBarPlot <- function(file_path_1, sheet, is.num = FALSE) {

	dat1 <- read.xlsx(xlsxFile = file_path_1, sheet = sheet, skipEmptyRows = TRUE)

	colnames(dat1)[1] <- "Experiment"
	colnames(dat1)[2] <- "Cluster"

	if (is.num) {
		dat1$Cluster <- factor(dat1$Cluster, levels = sort(unique(as.numeric(as.character(dat1$Cluster)))))
	}

	dat1 <- dat1 %>%
        group_by(Experiment) %>%
        mutate(Percentage = round(100*Cell_Count/sum(Cell_Count, na.rm = TRUE), digits=2)) %>%
        as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

	p1 <- ggplot(dat1, aes(fill=Experiment, y=Cell_Count, x=Cluster)) +
		geom_bar(position="stack", stat="identity") +
		geom_text(aes(label=Cell_Count), position=position_stack(vjust=0.5), size=3) +
		labs(title=paste0(sheet)) +
		theme(axis.text.x = element_text(angle=45, hjust=1))

	p2 <- ggplot(dat1, aes(fill=Experiment, y=Percentage, x=Cluster)) +
        geom_bar(position="stack", stat="identity") +
        geom_text(aes(label=Percentage), position=position_stack(vjust=0.5), size=3) +
        labs(title=paste0(sheet)) +
        theme(axis.text.x = element_text(angle=45, hjust=1))

    p <- p1 + p2

	ggsave(
		filename = paste0("cell_counts_", sheet, ".png"),
		plot = p,
		path = output_path,
		width = 21,
		height = 7
	)
}

plotStackBarPlot(dat_file_path, "seurat_clusters", is.num = TRUE)
plotStackBarPlot(dat_file_path, "HumanPrimaryCellAtlasData")
plotStackBarPlot(dat_file_path, "BlueprintEncodeData")
plotStackBarPlot(dat_file_path, "MouseRNAseqData")
plotStackBarPlot(dat_file_path, "ImmGenData")
plotStackBarPlot(dat_file_path, "DbImmuneCellExpressionData")
plotStackBarPlot(dat_file_path, "NovershternHematopoieticData")
plotStackBarPlot(dat_file_path, "MonacoImmuneData")
