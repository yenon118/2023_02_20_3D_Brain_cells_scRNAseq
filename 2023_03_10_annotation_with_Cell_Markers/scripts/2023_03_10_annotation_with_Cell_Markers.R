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

folder_path = file.path("../../2023_03_06_extract_normalized_counts/output")

processed_dat <- readRDS(file = file.path(folder_path, paste0("data_", condition, ".rds")))


folder_path <- file.path("../../2023_03_03_annotate_cell_types/data")

dat_file_path <- file.path(folder_path, "Cell_marker_Human.xlsx")

dat_sheet_names <- getSheetNames(dat_file_path)
print(dat_sheet_names)

ref <- read.xlsx(xlsxFile = dat_file_path, sheet = dat_sheet_names[1], skipEmptyRows = TRUE)


##################################################
# Process data
##################################################

dat <- processed_dat %>%
    pivot_longer(!Gene, names_to = "Cell", values_to = "Measurement") %>%
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


df <- ref %>%
    select(cell_name, marker) %>%
    group_by(cell_name) %>%
    mutate(Count = n()) %>%
    ungroup() %>%
    filter(Count > 18) %>%
    select(cell_name, marker) %>%
    left_join(dat, by = c("marker" = "Gene")) %>%
    drop_na(Cell) %>%
    drop_na(Measurement) %>%
    mutate(Flag = ifelse(Measurement > 0, TRUE, FALSE)) %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

rate_df <- df %>%
    group_by(cell_name, Cell) %>%
    summarize(Rate = sum(Flag, na.rm = TRUE)/n()) %>%
    ungroup() %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

max_rate_df <- rate_df %>%
    group_by(Cell) %>%
    slice_max(Rate) %>%
    ungroup() %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

summ_dat <- max_rate_df %>%
    group_by(Cell) %>%
    summarize(Count = n()) %>%
    ungroup() %>%
    arrange(desc(Count)) %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

annotation_df <- max_rate_df
annotation_df$cell_name[annotation_df$Cell %in% summ_dat$Cell[summ_dat$Count > 1]] <- "unknown"
annotation_df <- annotation_df %>%
    select(cell_name, Cell) %>%
    distinct(.keep_all = TRUE) %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

full_cells <- colnames(processed_dat)[!str_detect(colnames(processed_dat), fixed("gene", ignore_case = TRUE))]
partial_cells <- full_cells[!(full_cells %in% annotation_df$Cell)]

if (length(partial_cells) > 0) {
    annotation_df <- rbind(
        annotation_df,
        data.frame(
            cell_name = "unknown",
            Cell = partial_cells,
            check.names = FALSE,
            stringsAsFactors = FALSE
        )
    )
}

print(head(df))
print(tail(df))
print(dim(df))

print(head(rate_df))
print(tail(rate_df))
print(dim(rate_df))

print(head(annotation_df))
print(tail(annotation_df))
print(dim(annotation_df))


##################################################
# Save data
##################################################

write.table(
    x = annotation_df,
    sep = "\t",
    file = file.path(output_path, paste0("annotation_", condition, ".txt")),
    na = "",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)
