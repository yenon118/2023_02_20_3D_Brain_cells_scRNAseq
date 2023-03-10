#!/usr/bin/Rscript --vanilla
rm(list=ls())

set.seed(1)

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(argparse)

library(openxlsx)

library(SCINA)


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

processed_dat <- processed_dat %>%
    column_to_rownames(var = "Gene") %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

ref <- ref %>%
    filter(str_detect(tissue_class, fixed("Brain", ignore_case=FALSE))) %>%
    filter(str_detect(cancer_type, regex("normal", ignore_case=TRUE))) %>%
    filter(str_detect(cell_type, regex("normal", ignore_case=TRUE))) %>%
    filter(!str_detect(cell_name, fixed("Lake et al.Science", ignore_case=TRUE))) %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)


signatures <- list()

unique_cell_type <- sort(unique(sort(ref$cell_name)))

for (i in 1:length(unique_cell_type)) {
	marker_list <- sort(unique(sort(c(
        ref$marker[ref$cell_name == unique_cell_type[i]],
        ref$Symbol[ref$cell_name == unique_cell_type[i]]
    ))))

    marker_list <- marker_list[!is.na(marker_list)]
    marker_list <- marker_list[marker_list != ""]

    if(length(marker_list) > 18) {
        signatures[[unique_cell_type[i]]] <- marker_list
    }
}

# print(head(processed_dat))
# print(tail(processed_dat))
# print(dim(processed_dat))

print(head(ref))
print(tail(ref))
print(dim(ref))


print(signatures)


scina_object <- SCINA(
	processed_dat,
	signatures, 
	max_iter = 500, 
	convergence_n = 10, 
	convergence_rate = 0.999, 
	sensitivity_cutoff = 0.9, 
	rm_overlap=TRUE, 
	allow_unknown=TRUE
)

print(str(scina_object))

print(head(scina_object$cell_labels))


dat <- data.frame(
    Cell = colnames(processed_dat),
    Annotation = scina_object$cell_labels,
    check.names = FALSE,
    stringsAsFactors = FALSE
)


##################################################
# Save data
##################################################

write.table(
    x = dat,
    sep = "\t",
    file = file.path(output_path, paste0("scina_annotation_", condition, ".txt")),
    na = "",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)
