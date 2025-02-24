rm(list=ls())
# load libraries
library(Seurat)
library(SeuratObject)
library(dplyr)
library(data.table)
library(readr)
library(tidyr)
    

## Load ATAC data
# Input path for ATAC files
input_path <- "/mnt/DATA3/timo/ATAC/atac_score_matrices_all_cell_types"
files <- list.files(input_path, full.names = TRUE)

# Looop through files and add them to a list
list_of_dfs <- lapply(files, function(file) {
  read.table(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
})

# Merge all dfs in the list
atac <- do.call(rbind, list_of_dfs)
rownames(atac) <- atac$GeneID

# remove GeneID column
atac <- atac %>% select(-GeneID)

## Load RNA data
rna <- read.table("/mnt/DATA3/timo/data/sngl_cell_Tabula_Sapiens/avg_expression_Tabula_Sapiens_cell_type.tsv", 
                                   sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# Retain common ENSG, that are in atac
common_ensg <- intersect(rownames(atac), rownames(rna))
length(common_ensg)
# Subset both data frames to retain only common row names
rna_filtered <- rna[common_ensg, , drop = FALSE]
atac_filtered <- atac[common_ensg, , drop = FALSE]

# Remove 'RNA.' from the column names
colnames(rna_filtered) <- gsub("^RNA\\.", "", colnames(rna_filtered))

# Replace '.' with a space in the column names
colnames(rna_filtered) <- gsub("\\.", " ", colnames(rna_filtered))


atac_final <- atac_filtered
rna_final <- rna_filtered
    
print("Save outputfile.")
write.table(rna_final,"rna_all_cell_types.tsv", sep = "\t", row.names = TRUE, quote = FALSE)
write.table(atac_final, "atac_all_cell_types.tsv", sep = "\t", row.names = TRUE, quote = FALSE)


