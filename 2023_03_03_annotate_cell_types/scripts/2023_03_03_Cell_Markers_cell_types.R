#!/usr/bin/Rscript --vanilla
rm(list=ls())

set.seed(1)

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(argparse)

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
# Input and reference files
##################################################

folder_path = file.path(paste0("../../2023_02_23_Seurat_", condition, "/output/SignificantMarkers"))

dat <- read.table(
    file = file.path(folder_path, "data_markers_seurat_clusters.txt"),
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = ""
)


folder_path <- file.path("../data")

dat_file_path <- file.path(folder_path, "Cell_marker_Human.xlsx")

dat_sheet_names <- getSheetNames(dat_file_path)
print(dat_sheet_names)

ref <- read.xlsx(xlsxFile = dat_file_path, sheet = dat_sheet_names[1], skipEmptyRows = TRUE)


##################################################
# Process data
##################################################

dat$marker <- dat$gene

dat <- dat %>%
    group_by(cluster) %>%
    mutate(Total_Count = n_distinct(gene)) %>%
    ungroup() %>%
    mutate(cluster = paste0(cluster, " (", Total_Count, ")")) %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

ref <- ref %>%
    filter(str_detect(tissue_class, fixed("Brain", ignore_case=FALSE))) %>%
    filter(str_detect(cancer_type, regex("normal", ignore_case=TRUE))) %>%
    filter(str_detect(cell_type, regex("normal", ignore_case=TRUE))) %>%
    filter(!str_detect(cell_name, fixed("Lake et al.Science", ignore_case=TRUE))) %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)


print(head(dat))
print(tail(dat))
print(dim(dat))

print(head(ref))
print(tail(ref))
print(dim(ref))


df <- dat %>%
    left_join(ref, by = "marker") %>%
    drop_na(cell_name) %>%
    filter(cell_name != "") %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)


print(head(df))
print(tail(df))
print(dim(df))


##################################################
# Plotting
##################################################

df_plotting <- df %>%
    group_by(cell_name, cluster) %>%
    summarize(Count = n_distinct(marker)) %>%
    ungroup() %>%
    arrange(cluster, Count, cell_name) %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

print(head(df_plotting))
print(tail(df_plotting))
print(dim(df_plotting))


p <- ggplot(df_plotting, aes(x=cell_name, y=Count, fill=cell_name)) +
    geom_bar(stat='identity') +
    geom_text(aes(label = Count), vjust=-0.3, stat="identity", size=3) +
    facet_wrap(~cluster, scales="free_y", ncol=3) +
    theme(
        axis.text.x = element_text(angle=75, hjust=1), 
        legend.position="none"
    )


ggsave(
    filename = paste0("Cell_Markers_cell_types_", condition, ".png"),
    plot = p,
    path = output_path,
    width = 16,
    height = 12,
    units = "in",
    dpi = 300
)

