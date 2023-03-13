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

dat_file_path <- file.path(folder_path, "PanglaoDB_markers_27_Mar_2020.tsv")

ref <- read.table(
    file = dat_file_path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
)


##################################################
# Process data
##################################################

dat$official.gene.symbol <- dat$gene

dat <- dat %>%
    group_by(cluster) %>%
    mutate(Total_Count = n_distinct(gene)) %>%
    ungroup() %>%
    mutate(cluster = paste0(cluster, " (", Total_Count, ")")) %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

ref <- ref %>%
    filter(str_detect(species, fixed("Hs", ignore_case=FALSE))) %>%
    filter(str_detect(organ, fixed("brain", ignore_case=TRUE))) %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

temp_ref <- ref 
temp_ref$official.gene.symbol <- temp_ref$nicknames

temp_ref <- temp_ref %>%
    separate_rows(official.gene.symbol, sep = "\\|", convert = TRUE) %>%
    distinct(.keep_all = TRUE) %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

ref <- rbind(ref, temp_ref)

ref <- ref %>%
    distinct(.keep_all = TRUE) %>%
    distinct(official.gene.symbol, .keep_all = TRUE) %>%
    filter(cell.type != "") %>%
    filter(official.gene.symbol != "") %>%
    drop_na(cell.type) %>%
    drop_na(official.gene.symbol) %>%
    distinct(.keep_all = TRUE) %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)


print(head(dat))
print(tail(dat))
print(dim(dat))

print(head(ref))
print(tail(ref))
print(dim(ref))


df <- dat %>%
    left_join(ref, by = "official.gene.symbol") %>%
    distinct(.keep_all = TRUE) %>%
    drop_na(cell.type) %>%
    filter(cell.type != "") %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)


print(head(df))
print(tail(df))
print(dim(df))


##################################################
# Save data
##################################################

write.table(
    x = df,
    file = file.path(output_path, paste0("PanglaoDB_cell_types_", condition, ".txt")),
    sep = "\t",
    na = "",
    quote = FALSE,
    row.names = FALSE
)


##################################################
# Plotting
##################################################

df_plotting <- df %>%
    group_by(cell.type, cluster) %>%
    summarize(Count = n_distinct(official.gene.symbol, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(cluster, Count, cell.type) %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

print(head(df_plotting))
print(tail(df_plotting))
print(dim(df_plotting))


p <- ggplot(df_plotting, aes(x=cell.type, y=Count, fill=cell.type)) +
    geom_bar(stat='identity') +
    geom_text(aes(label = Count), vjust=-0.3, stat="identity", size=3) +
    facet_wrap(~cluster, scales="free_y", ncol=3) +
    theme(
        axis.text.x = element_text(angle=75, hjust=1), 
        legend.position="none"
    )


ggsave(
    filename = paste0("PanglaoDB_cell_types_", condition, ".png"),
    plot = p,
    path = output_path,
    width = 16,
    height = 12,
    units = "in",
    dpi = 300
)
