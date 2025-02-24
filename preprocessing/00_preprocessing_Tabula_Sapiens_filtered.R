#sessionInfo()  # See all loaded packages and their versions
library(Matrix)
library(Seurat)
library(SeuratObject)

# Load data
tabula_sapiens <- readRDS("/mnt/DATA3/timo/data/sngl_cell_Tabula_Sapiens/981bcf57-30cb-4a85-b905-e04373432fef.rds")

dim <- dim(tabula_sapiens)
print(paste("Initial Dimension: ", dim[1], dim[2]))

# Subset for assay = 10x 3' v3
subset_tabula_sapiens <- subset(tabula_sapiens, subset = assay == "10x 3' v3")

dim_subset <- dim(subset_tabula_sapiens)
print(paste("Filtered for assay = 10x 3' v3, Dimension:", dim_subset[1], dim_subset[2]))

# Subset the data by excluding cells with "sperm" in the "cell type" column
subset_tabula_sapiens <- subset(subset_tabula_sapiens, subset = cell_type != "sperm")
dim_subset <- dim(subset_tabula_sapiens)
print(paste("Filtered cell_type != sperm, Dimension: ", dim_subset[1], dim_subset[2]))

# Create additional column "cell_type_tissue" by concatenating "cell_type" and "tissue_in_publication"
subset_tabula_sapiens$cell_type_tissue <- paste(subset_tabula_sapiens$cell_type, subset_tabula_sapiens$tissue_in_publication, sep = "_")

dim_final <- dim(tabula_sapiens)
print(paste("Final Dimension: ", dim_final[1], dim_final[2]))

# Save filtered Seurat Object
saveRDS(subset_tabula_sapiens, file = "/mnt/DATA3/timo/data/sngl_cell_Tabula_Sapiens/Tabula_Sapiens_filtered.rds")