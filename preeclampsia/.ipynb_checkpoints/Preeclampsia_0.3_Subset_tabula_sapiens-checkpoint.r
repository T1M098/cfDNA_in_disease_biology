#https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

#install.packages("devtools")
#devtools::install_github("immunogenomics/presto")

#remove.packages("Matrix")
#library(devtools)
#install_version("Matrix", version = "1.6.4")

rm(list=ls())
setwd("/mnt/DATA3/timo")
getwd()
# https://github.com/immunogenomics/presto
# Find markers using presto
#install.packages("devtools")
#devtools::install_github("immunogenomics/presto")
# load libraries
library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(presto)
library(future)
library(data.table)
library(ggplot2)
#plan(multicore, workers = 10)
options(future.globals.maxSize = 200*1024^3)

sessionInfo()

packageVersion("Seurat")
packageVersion("Matrix")

print("Loading .rds test file.")
tabula_sapiens <- readRDS("/mnt/DATA3/timo/data/sngl_cell_Tabula_Sapiens/Tabula_Sapiens_filtered.rds") # Tabula_Sapiens_filtered.rds seurat_subset_10k_test.rds
print("File loaded.")

# remove unknown cell_type
tabula_sapiens <- subset(tabula_sapiens, subset = cell_type != "unknown")

tabula_sapiens@assays

tabula_sapiens <- FindVariableFeatures(tabula_sapiens, selection.method = "vst", nfeatures = 2000)

scale_in_chunks <- function(seurat_object, genes, chunk_size = 1000) {
  # Initialize an empty scale.data matrix if not already present
  if (sum(dim(tabula_sapiens@assays$RNA@scale.data) == c(0, 0)) == 2) {
    seurat_object@assays$RNA@scale.data <- matrix(NA, nrow = length(genes), ncol = ncol(seurat_object),
                                                  dimnames = list(genes, colnames(seurat_object)))
  }
  
  # Loop through the genes in chunks
  for (i in seq(1, length(genes), by = chunk_size)) {
    chunk_genes <- genes[i:min(i + chunk_size - 1, length(genes))]
    
    # Scale the chunk of genes
    temp_object <- ScaleData(seurat_object, features = chunk_genes, verbose = FALSE)
    
    # Extract the scaled data for the chunk and ensure it has the correct structure
    scaled_chunk <- temp_object@assays$RNA@scale.data[chunk_genes, , drop = FALSE]
    
    # Update the scale.data slot incrementally
    seurat_object@assays$RNA@scale.data[chunk_genes, ] <- scaled_chunk
  }
  
  return(seurat_object)
}

DefaultAssay(tabula_sapiens) <- "RNA"

all.genes <- rownames(tabula_sapiens)

tabula_sapiens <- ScaleData(tabula_sapiens, features = VariableFeatures(tabula_sapiens), verbose = FALSE)

#tabula_sapiens <- scale_in_chunks(tabula_sapiens, all.genes)

tabula_sapiens <- RunPCA(tabula_sapiens, features = VariableFeatures(tabula_sapiens), dims = 1:20)

# Examine and visualize PCA results a few different ways
print(tabula_sapiens[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(tabula_sapiens, reduction = "pca", raster=FALSE) + NoLegend()

tabula_sapiens <- FindNeighbors(tabula_sapiens, dims = 1:10)
tabula_sapiens <- FindClusters(tabula_sapiens, resolution = 0.1)

tabula_sapiens <- RunUMAP(tabula_sapiens, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(tabula_sapiens, reduction = "umap", label = TRUE) + ggtitle("UMAP of original Tabula Sapiens Clusters")

#saveRDS(tabula_sapiens, file = "/mnt/DATA3/timo/data/sngl_cell_Tabula_Sapiens/Tabula_Sapiens_original_clusters.rds")

tabula_sapiens <- readRDS(file = "/mnt/DATA3/timo/data/sngl_cell_Tabula_Sapiens/Tabula_Sapiens_original_clusters.rds")

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(tabula_sapiens, reduction = "umap", label = TRUE)

## remove clusters
# 24 : myoepithelial cell
# 29: eye photoreceptors
tabula_sapiens <- subset(tabula_sapiens, subset = seurat_clusters != 29)
tabula_sapiens <- subset(tabula_sapiens, subset = seurat_clusters != 24)

# merge cluster 7, 16 and 21 to 7 ("exocrine cells")
Idents(tabula_sapiens)[Idents(tabula_sapiens) %in% c("7", "16", "21")] <- "7"

clusters_df_tabula <- data.frame(Cell = rownames(tabula_sapiens@meta.data),
                          cell_type = tabula_sapiens@meta.data$cell_type,
                          cell_type_tissue = tabula_sapiens@meta.data$cell_type_tissue,
                          Cluster = tabula_sapiens@meta.data$seurat_clusters)

# merge cluster 7,16 and 21 => exocrine cells
clusters_df_tabula <- clusters_df_tabula %>%
  mutate(Cluster = as.numeric(as.character(Cluster))) %>%
  mutate(Cluster_new = case_when(
    Cluster %in% c(16,21) ~ 7,      # Assign 7 when Cluster is 16 or 21
    TRUE ~ Cluster           # Keep the original Cluster value for other cases
  ))

new.cluster.ids_tabula <- list()
for (i in sort(unique(clusters_df_tabula$Cluster_new))){
    i <- as.numeric(i)
    cell_type_names <- as.character((clusters_df_tabula[clusters_df_tabula$Cluster_new == i,])$cell_type)
    cell_type_counts <- table(cell_type_names)
    #print(paste("Cluster: ", i))
    most_abundant_cell_type <- names(as.list(head(sort(cell_type_counts, decreasing = TRUE), 1)))
    top3 <- names(as.list(head(sort(cell_type_counts, decreasing = TRUE), 3)))
    #print(top3)
    new.cluster.ids_tabula[i+1] <- most_abundant_cell_type    
    }
cluster_cell_types <- unlist(new.cluster.ids_tabula)


# Create a mapping table
mapping_table <- data.frame(
  Cluster_new = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 22, 23, 25, 26, 27, 28),
  Cluster_cell_type = cluster_cell_types)


# Join with the main dataframe
clusters_df_tabula <- clusters_df_tabula %>%
  left_join(mapping_table, by = "Cluster_new")


# add new cell type names
# exocrine cells
clusters_df_tabula <- clusters_df_tabula %>%
  mutate(Cluster_new_cell_type = case_when(
    Cluster_new == 7 ~ "exocrine cells",  # If Cluster_new is 16, assign "exocrine cells"
    Cluster_new == 0 ~ "B, T and NK cells",
    TRUE ~ Cluster_cell_type  # Otherwise, retain the original Cluster_cell_type value
  ))

final.cluster.ids_tabula <- list()
for (i in sort(unique(clusters_df_tabula$Cluster_new))){
    cell_type <- clusters_df_tabula %>% filter(Cluster_new == i) %>% pull(Cluster_new_cell_type) %>% unique()
    # Add the cell type(s) to the list with the cluster number as the name
    final.cluster.ids_tabula[[as.character(i)]] <- cell_type
}

names(final.cluster.ids_tabula) <- levels(tabula_sapiens)
tabula_sapiens <- RenameIdents(tabula_sapiens, final.cluster.ids_tabula)

DimPlot(tabula_sapiens, reduction = "umap", label = TRUE, label.size = 3, pt.size = 0.5) + NoLegend() + ggtitle("UMAP of Tabula Sapiens New Clusters")

for (i in sort(unique(clusters_df_tabula$Cluster_new))){
    cell_types <- as.character((clusters_df_tabula[clusters_df_tabula$Cluster == i,])$cell_type)
    sorted_cell_types <- sort(table(cell_types), decreasing = TRUE)
    cluster_cell_type <- clusters_df_tabula %>% filter(Cluster == i) %>% pull(Cluster_new_cell_type) %>% unique()
    print(paste("Cluster: ", i, cluster_cell_type))
    print(sorted_cell_types[1:10])
    }

saveRDS(tabula_sapiens, file = "/mnt/DATA3/timo/data/sngl_cell_Tabula_Sapiens/Tabula_Sapiens_new_clusters.rds")

tabula_sapiens <- readRDS(file = "/mnt/DATA3/timo/data/sngl_cell_Tabula_Sapiens/Tabula_Sapiens_new_clusters.rds")

DimPlot(tabula_sapiens, reduction = "umap", label = TRUE, label.size = 3, pt.size = 0.5) + NoLegend() + ggtitle("UMAP of Tabula Sapiens New Clusters")

length(unique(Idents(tabula_sapiens)))

## subset and take 2000k cells from each cluster

sample_size <- 3000

head(tabula_sapiens@meta.data)

set.seed(42)
# Get the cluster identities
cluster_ids <- unique(Idents(tabula_sapiens))

# Initialize sampled_cells as a list
sampled_cells <- list()
s <- 0  # Total number of cells processed

for (cluster_id in cluster_ids) {
    # Subset cells in the current cluster
    cells_in_cluster <- WhichCells(tabula_sapiens, idents = cluster_id)
    
    # Sample or use all cells if fewer than sample_size
    if (length(cells_in_cluster) < sample_size) {
        sampled_cells[[as.character(cluster_id)]] <- cells_in_cluster
        print(cluster_id)
    } else {
        sampled_cells[[as.character(cluster_id)]] <- sample(cells_in_cluster, sample_size)
    }
    
    # Print length of cells in the cluster and accumulate total
    print(length(cells_in_cluster))
    s <- s + length(cells_in_cluster)
}

# Print summary
print(s)



# Convert the list of cell names to a vector
cells_to_keep <- unlist(sampled_cells)
str(cells_to_keep)

# Subset the Seurat object to keep only these cells
tabula_sapiens_subset <- subset(tabula_sapiens, cells = cells_to_keep)

tabula_sapiens_subset@assays

tabula_sapiens_subset <- FindVariableFeatures(tabula_sapiens_subset, selection.method = "vst", nfeatures = 2000)
DefaultAssay(tabula_sapiens_subset) <- "RNA"

tabula_sapiens_subset <- ScaleData(tabula_sapiens_subset, features = VariableFeatures(tabula_sapiens_subset), verbose = FALSE)
tabula_sapiens_subset <- RunPCA(tabula_sapiens_subset, features = VariableFeatures(tabula_sapiens_subset), dims = 1:20)
# Examine and visualize PCA results a few different ways
print(tabula_sapiens_subset[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(tabula_sapiens_subset, reduction = "pca", raster=FALSE) + NoLegend()
tabula_sapiens_subset <- FindNeighbors(tabula_sapiens_subset, dims = 1:10)
tabula_sapiens_subset <- FindClusters(tabula_sapiens_subset, resolution = 0.1)
tabula_sapiens_subset <- RunUMAP(tabula_sapiens_subset, dims = 1:10)

DimPlot(tabula_sapiens_subset, reduction = "umap", label = TRUE, label.size = 3, pt.size = 0.5) + NoLegend() + ggtitle("UMAP of 3k Subset Tabula Sapiens New Clusters")

clusters_df_tabula_2 <- data.frame(Cell = rownames(tabula_sapiens_subset@meta.data),
                          cell_type = tabula_sapiens_subset@meta.data$cell_type,
                          cell_type_tissue = tabula_sapiens_subset@meta.data$cell_type_tissue,
                          Cluster = tabula_sapiens_subset@meta.data$seurat_clusters)

new.cluster.ids_tabula_sapiens_subset <- list()
for (i in 0:(length(unique(clusters_df_tabula_2$Cluster))-1)){
    #cell_type_names <- as.character((clusters_df[clusters_df$Cluster == i,])$cell_type_tissue)
    #cell_type_counts <- table(cell_type_names)
    #print(paste("Cluster: ", i))
    #print(head(sort(cell_type_counts, decreasing = TRUE), 3))
    cell_type_names <- as.character((clusters_df_tabula_2[clusters_df_tabula_2$Cluster == i,])$cell_type)
    cell_type_counts <- table(cell_type_names)
    print(paste("Cluster: ", i))
    most_abundant_cell_type <- names(as.list(head(sort(cell_type_counts, decreasing = TRUE), 1)))
    top3 <- names(as.list(head(sort(cell_type_counts, decreasing = TRUE), 3)))
    print(top3)
    new.cluster.ids_tabula_sapiens_subset <- c(new.cluster.ids_tabula_sapiens_subset, most_abundant_cell_type)
    #if (most_abundant_cell_type %in% new.cluster.ids_placenta){
    #    most_abundant_cell_type <- names(as.list(head(sort(cell_type_counts, decreasing = TRUE), 2)))[2]
    #    new.cluster.ids_placenta <- c(new.cluster.ids_placenta, most_abundant_cell_type)
    #} else{
    #    new.cluster.ids_placenta <- c(new.cluster.ids_placenta, most_abundant_cell_type)
    #}
    
    }

new.cluster.ids_tabula_sapiens_subset

names(new.cluster.ids_tabula_sapiens_subset) <- levels(tabula_sapiens_subset)
tabula_sapiens_subset <- RenameIdents(tabula_sapiens_subset, new.cluster.ids_tabula_sapiens_subset)

DimPlot(tabula_sapiens_subset, reduction = "umap", label = TRUE, label.size = 3, pt.size = 0.5) + NoLegend() + ggtitle("UMAP of 3k Subset Tabula Sapiens original Clusters") + theme(plot.title = element_text(size = 12))

saveRDS(tabula_sapiens_subset, file = "/mnt/DATA3/timo/data/sngl_cell_Tabula_Sapiens/Tabula_Sapiens_subset3k_original_clusters.rds")

tabula_sapiens_subset <- readRDS("/mnt/DATA3/timo/data/sngl_cell_Tabula_Sapiens/Tabula_Sapiens_subset3k_original_clusters.rds")

levels(tabula_sapiens_subset)

head(Idents(tabula_sapiens_subset))

# library(ggplot2)
umap_plot <- DimPlot(tabula_sapiens_subset, reduction = "umap", label = TRUE, label.size = 4, pt.size = 0.5) +
    NoLegend() + ggtitle("UMAP Clustering of Tabula Sapiens Subset") +
    theme(plot.title = element_text(size = 12)) +
    labs(x = "UMAP 1", y = "UMAP 2")
umap_plot

ggsave("/mnt/DATA3/timo/plots/umap_plot_tabula_sapiens_subset.pdf", plot = umap_plot, width = 10, height = 8)
