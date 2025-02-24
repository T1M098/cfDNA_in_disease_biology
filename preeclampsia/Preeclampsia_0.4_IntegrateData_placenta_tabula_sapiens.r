rm(list=ls())
setwd("/mnt/DATA3/timo")
getwd()
# increase memory usage

library(Matrix)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(Seurat)
library(SeuratObject)
library(data.table)
library(dplyr)

tabula_sapiens <- readRDS("/mnt/DATA3/timo/data/sngl_cell_Tabula_Sapiens/Tabula_Sapiens_subset3k_original_clusters.rds")
print("tabula_sapiens loaded.")
placenta <- readRDS("/mnt/DATA3/timo/data/placenta/placenta_subset3k_original_clusters.rds")
print("placenta loaded.")
fetal_CA <- readRDS("/mnt/DATA3/timo/data/placenta/fetal_CA_subset3k_merged_clusters.rds")
print("fetal_CA loaded.")

tabula_sapiens <- UpdateSeuratObject(tabula_sapiens)

tabula_sapiens@assays
placenta@assays
fetal_CA@assays

colnames(tabula_sapiens@meta.data)

colnames(placenta@meta.data)

colnames(fetal_CA@meta.data)

# testing
#tabula_sapiens <- subset(tabula_sapiens, cells = sample(Cells(tabula_sapiens), 1000))
#placenta <- subset(placenta, cells = sample(Cells(placenta), 1000))
#fetal_CA <- subset(fetal_CA, cells = sample(Cells(fetal_CA), 1000))

print("tabula_sapiens FindVarF.")
tabula_sapiens <- FindVariableFeatures(tabula_sapiens, selection.method = "vst", nfeatures = 2000)
print("placenta FindVarF.")
placenta <- FindVariableFeatures(placenta, selection.method = "vst", nfeatures = 2000)
print("fetal_CA FindVarF.")
fetal_CA <- FindVariableFeatures(fetal_CA, selection.method = "vst", nfeatures = 2000)

DefaultAssay(tabula_sapiens) <- "RNA"
DefaultAssay(placenta) <- "RNA"
DefaultAssay(fetal_CA) <- "RNA"

length(VariableFeatures(tabula_sapiens))
length(VariableFeatures(placenta))
length(VariableFeatures(fetal_CA))

anchors <- FindIntegrationAnchors(object.list = list(tabula_sapiens, placenta, fetal_CA), reference = 1, dim = 1:10)
#anchors <- FindIntegrationAnchors(object.list = list(tabula_sapiens, fetal_cell_atlas_placenta_seurat_object, placenta_seurat_object), reference = 1, dims = 1:10) #

print("Anchors: ")
print(anchors)

# combine them
print("Integrating Datasets into Tabula Sapiens")
integrated_seurat_fetal_and_placenta <- IntegrateData(anchorset = anchors)

integrated_seurat_fetal_and_placenta@assays

unique(integrated_seurat_fetal_and_placenta@meta.data[, c("tissue")])

# Replace NA values in 'tissue_in_publication' with 'tissue' => either decidua or placenta
integrated_seurat_fetal_and_placenta@meta.data <- integrated_seurat_fetal_and_placenta@meta.data %>%
  mutate(tissue_in_publication = ifelse(is.na(tissue_in_publication), tissue, tissue_in_publication))

# Redo cell_type_tissue column
integrated_seurat_fetal_and_placenta@meta.data$cell_type_tissue <- paste(
  integrated_seurat_fetal_and_placenta@meta.data$cell_type, 
  integrated_seurat_fetal_and_placenta@meta.data$tissue_in_publication, 
  sep = "_"
)


saveRDS(integrated_seurat_fetal_and_placenta, file = "/mnt/DATA3/timo/data/placenta/Integrated_data.rds")

tabula_sapiens_placenta <- readRDS("/mnt/DATA3/timo/data/placenta/Integrated_data.rds")
print("File loaded.")
print(tabula_sapiens_placenta@assays)

##################################
print("Average Expression per Cell Type Tissue")
avg_expr <- AverageExpression(object = tabula_sapiens_placenta,
                       assays = "RNA",
                       group.by = "cell_type_tissue",
                       return.seurat = FALSE)
avg_expr_df <- as.data.frame(avg_expr)
print("Only keep Genes with >= 3 non-zero values.")
# Count non-zero entries for each row directly on the sparse matrix
non_zero_counts <- rowSums(avg_expr_df != 0)
# Create a logical vector to keep rows with at least 3 or more non-zero values
threshold <- 3
keep_rows <- non_zero_counts[non_zero_counts > threshold]
rows_names_keep <- names(keep_rows)
# subset the initial dgCMatrix using the row names (ENSG00010123) which contain at least 3 non-zero values
avg_expr_df_filtered <- avg_expr_df[rownames(avg_expr_df) %in% rows_names_keep, ]



# save file
write.table(avg_expr_df_filtered, file = "/mnt/DATA3/timo/data/sngl_cell_Tabula_Sapiens/avg_expression_Tabula_Sapiens_cell_type_tissue_placenta_new.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
