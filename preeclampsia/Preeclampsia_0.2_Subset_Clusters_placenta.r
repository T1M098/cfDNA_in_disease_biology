#https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

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
library(patchwork)
library(data.table)
library(presto)
library(ggplot2)
library(future)
#plan(multicore, workers = 10)
options(future.globals.maxSize = 250*1024^3)

packageVersion("Seurat")

# Define a function to scale in chunks
scale_in_chunks <- function(seurat_object, genes, chunk_size = 1000) {
  for (i in seq(1, length(genes), by = chunk_size)) {
    chunk_genes <- genes[i:min(i + chunk_size - 1, length(genes))]
    seurat_object <- ScaleData(seurat_object, features = chunk_genes, verbose = FALSE)
  }
  return(seurat_object)
}

print("Loading .rds test file.")
placenta <- readRDS("/mnt/DATA3/timo/data/placenta/placenta_seurat_object_norm.rds") # Tabula_Sapiens_filtered.rds seurat_subset_10k_test.rds
print("File loaded.")

head(placenta@meta.data)

placenta@assays

placenta <- FindVariableFeatures(placenta, selection.method = "vst", nfeatures = 2000)

DefaultAssay(placenta) <- "RNA"

#all.genes <- rownames(placenta)
#placenta <- scale_in_chunks(placenta, all.genes)

placenta <- ScaleData(placenta, features = rownames(placenta))

#placenta <- RunPCA(placenta)
placenta <- RunPCA(placenta, features = VariableFeatures(placenta), dims = 1:20)

# Examine and visualize PCA results a few different ways
print(placenta[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(placenta, reduction = "pca") + NoLegend()

placenta <- FindNeighbors(placenta, dims = 1:10)
placenta <- FindClusters(placenta, resolution = 0.1)

length(unique(placenta@meta.data$cell_type))

unique(placenta@meta.data$cell_type)

placenta <- RunUMAP(placenta, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(placenta, reduction = "umap", label = TRUE) + ggtitle("UMAP of original Placenta Clusters")

#saveRDS(placenta, file = "/mnt/DATA3/timo/data/placenta/placenta_original_clusters.rds")

placenta <- readRDS(file = "/mnt/DATA3/timo/data/placenta/placenta_original_clusters.rds")

# individual clusters
DimPlot(placenta, reduction = "umap", label = TRUE) + ggtitle("UMAP of original Placenta Clusters")

head(placenta@meta.data)

placenta@meta.data$cell_type_tissue <- paste(placenta@meta.data$cell_type, placenta@meta.data$tissue, sep = "_")

# EVT, SCT, VCT, fFB1, fFB2, HB
# extravillous cytotrophoblast (EVT), proliferative extravillous cytotrophoblasts (EVTp), syncytiotrophoblast (SCT), 
# proliferative villous cytotrophoblasts (VCTp),  villous cytotrophoblasts (VCT), Hoffbauer cells (HB), 2x fetal fibroblasts (fFB1 & 2)

unique(placenta@meta.data$cell_type)

clusters_df_placenta <- data.frame(Cell = rownames(placenta@meta.data),
                          cell_type = placenta@meta.data$cell_type,
                          cell_type_tissue = placenta@meta.data$cell_type_tissue,
                          Cluster = placenta$seurat_clusters)

cell_types <- as.character((clusters_df_placenta[clusters_df_placenta$Cluster == 5,])$cell_type_tissue)
sorted_cell_types <- sort(table(cell_types), decreasing = TRUE)
sorted_cell_types
cell_types <- as.character((clusters_df_placenta[clusters_df_placenta$Cluster == 2,])$cell_type_tissue)
sorted_cell_types <- sort(table(cell_types), decreasing = TRUE)
sorted_cell_types
cell_types <- as.character((clusters_df_placenta[clusters_df_placenta$Cluster == 7,])$cell_type_tissue)
sorted_cell_types <- sort(table(cell_types), decreasing = TRUE)
sorted_cell_types

new.cluster.ids_placenta <- c()
for (i in 0:(length(unique(clusters_df_placenta$Cluster))-1)){
    #cell_type_names <- as.character((clusters_df[clusters_df$Cluster == i,])$cell_type_tissue)
    #cell_type_counts <- table(cell_type_names)
    #print(paste("Cluster: ", i))
    #print(head(sort(cell_type_counts, decreasing = TRUE), 3))
    cell_type_names <- as.character((clusters_df_placenta[clusters_df_placenta$Cluster == i,])$cell_type)
    cell_type_counts <- table(cell_type_names)
    print(paste("Cluster: ", i))
    most_abundant_cell_type <- names(as.list(head(sort(cell_type_counts, decreasing = TRUE), 1)))
    top3 <- names(as.list(head(sort(cell_type_counts, decreasing = TRUE), 3)))
    print(top3)
    new.cluster.ids_placenta <- c(new.cluster.ids_placenta, most_abundant_cell_type)
    #if (most_abundant_cell_type %in% new.cluster.ids_placenta){
    #    most_abundant_cell_type <- names(as.list(head(sort(cell_type_counts, decreasing = TRUE), 2)))[2]
    #    new.cluster.ids_placenta <- c(new.cluster.ids_placenta, most_abundant_cell_type)
    #} else{
    #    new.cluster.ids_placenta <- c(new.cluster.ids_placenta, most_abundant_cell_type)
    #}
    
    }

new.cluster.ids_placenta <- c()
for (i in 0:(length(unique(clusters_df_placenta$Cluster))-1)){
    #cell_type_names <- as.character((clusters_df[clusters_df$Cluster == i,])$cell_type_tissue)
    #cell_type_counts <- table(cell_type_names)
    #print(paste("Cluster: ", i))
    #print(head(sort(cell_type_counts, decreasing = TRUE), 3))
    cell_type_names <- as.character((clusters_df_placenta[clusters_df_placenta$Cluster == i,])$cell_type_tissue)
    cell_type_counts <- table(cell_type_names)
    print(paste("Cluster: ", i))
    most_abundant_cell_type <- names(as.list(head(sort(cell_type_counts, decreasing = TRUE), 1)))
    top3 <- names(as.list(head(sort(cell_type_counts, decreasing = TRUE), 3)))
    print(top3)
    new.cluster.ids_placenta <- c(new.cluster.ids_placenta, most_abundant_cell_type)
    #if (most_abundant_cell_type %in% new.cluster.ids_placenta){
    #    most_abundant_cell_type <- names(as.list(head(sort(cell_type_counts, decreasing = TRUE), 2)))[2]
    #    new.cluster.ids_placenta <- c(new.cluster.ids_placenta, most_abundant_cell_type)
    #} else if (most_abundant_cell_type %in% new.cluster.ids_placenta == FALSE){
    #    new.cluster.ids_placenta <- c(new.cluster.ids_placenta, most_abundant_cell_type)
    #}
    
    }

unique(new.cluster.ids_placenta)

levels(placenta)

names(new.cluster.ids_placenta) <- levels(placenta)
placenta <- RenameIdents(placenta, new.cluster.ids_placenta)

DimPlot(placenta, reduction = "umap", label = TRUE, label.size = 3, pt.size = 0.5) + NoLegend() + ggtitle("UMAP of original Placenta Clusters")

#### sample 2k cells
placenta <- readRDS(file = "/mnt/DATA3/timo/data/placenta/placenta_original_clusters.rds")

cells_in_cluster <- WhichCells(placenta, idents = 11)
n_cells <- as.integer(length(cells_in_cluster))
print(paste(cluster_id, n_cells))
n_cells < 2000

sampled_cells <- list()
sample_size <- 3000
for (cluster_id in cluster_ids) {
    cells_in_cluster <- WhichCells(placenta, idents = cluster_id)
    n_cells <- as.integer(length(cells_in_cluster))
    if (n_cells < sample_size){
        sampled_cells[[as.character(cluster_id)]] <- cells_in_cluster
    } else {
        sampled_cells[[as.character(cluster_id)]] <- sample(cells_in_cluster, sample_size)
    }
    
}
str(sampled_cells)

# Convert the list of cell names to a vector
cells_to_keep <- unlist(sampled_cells)

str(cells_to_keep)

# Subset the Seurat object to keep only these cells
placenta_subset <- subset(placenta, cells = cells_to_keep)

placenta_subset <- FindVariableFeatures(placenta_subset, selection.method = "vst", nfeatures = 2000)

DefaultAssay(placenta_subset) <- "RNA"

placenta_subset <- ScaleData(placenta_subset, features = rownames(placenta))

placenta_subset <- RunPCA(placenta_subset, features = VariableFeatures(placenta), dims = 1:20)

DimPlot(placenta_subset, reduction = "pca") + NoLegend()

placenta_subset <- FindNeighbors(placenta_subset, dims = 1:10)
placenta_subset <- FindClusters(placenta_subset, resolution = 0.1)

placenta_subset <- RunUMAP(placenta_subset, dims = 1:10)

DimPlot(placenta_subset, reduction = "umap", label = TRUE) + ggtitle("UMAP of subset Placenta Clusters")

clusters_df_placenta_subset <- data.frame(Cell = rownames(placenta_subset@meta.data),
                          cell_type = placenta_subset@meta.data$cell_type,
                          cell_type_tissue = placenta_subset@meta.data$tissue,
                          Cluster = placenta_subset$seurat_clusters)

# Initialize the vector to store new cluster IDs
new.cluster.ids_placenta_subset_subset <- c()
used_cell_types <- c()  # Keep track of already assigned cell types

# Loop through each cluster
for (i in levels(clusters_df_placenta_subset$Cluster)){
    i <- as.numeric(i)
    
    # Get the cell types for the current cluster
    cell_type_names <- as.character(clusters_df_placenta_subset[clusters_df_placenta_subset$Cluster == i, ]$cell_type)
    
    # Count the occurrences of each cell type in the cluster
    cell_type_counts <- table(cell_type_names)
    
    # Get the top 3 most abundant cell types
    top3 <- names(sort(cell_type_counts, decreasing = TRUE)[1:3])
    # print(paste("Cluster: ", i))
    # print("Top 3 cell types:")
    # print(top3)
    
    # Check if the most abundant cell type is already used as a cluster ID
    chosen_cell_type <- top3[1]
    
    # Loop through the top 3 and select the first available one that hasn't been used yet
    for (cell_type in top3) {
        if (!(cell_type %in% used_cell_types)) {
            chosen_cell_type <- cell_type
            used_cell_types <- c(used_cell_types, chosen_cell_type)  # Mark this cell type as used
            break  # Exit the loop once a valid cell type is found
        }
    }
    
    # Append the chosen cell type to the new cluster ID list
    new.cluster.ids_placenta_subset_subset <- c(new.cluster.ids_placenta_subset_subset, chosen_cell_type)
}

new.cluster.ids_placenta_subset_subset


names(new.cluster.ids_placenta_subset_subset) <- levels(placenta_subset)
placenta_subset <- RenameIdents(placenta_subset, new.cluster.ids_placenta_subset_subset)

DimPlot(placenta_subset, reduction = "umap", label = TRUE, label.size = 4) + NoLegend() +  ggtitle("UMAP of Cell Clusters from Placenta Data Subset")

saveRDS(placenta_subset, file = "/mnt/DATA3/timo/data/placenta/placenta_subset3k_original_clusters.rds")

placenta_subset <- readRDS("/mnt/DATA3/timo/data/placenta/placenta_subset3k_original_clusters.rds")

DimPlot(placenta_subset, reduction = "umap", label = TRUE, label.size = 4) + ggtitle("UMAP Clustering Placenta Data Subset")

clusters_df_placenta_subset <- data.frame(Cell = rownames(placenta_subset@meta.data),
                          cell_type = placenta_subset@meta.data$cell_type,
                          cell_type_tissue = placenta_subset@meta.data$tissue,
                          Cluster = placenta_subset$seurat_clusters)

# Initialize the vector to store new cluster IDs
new.cluster.ids_placenta_subset_subset <- c()
used_cell_types <- c()  # Keep track of already assigned cell types

# Loop through each cluster
for (i in levels(clusters_df_placenta_subset$Cluster)){
    i <- as.numeric(i)
    
    # Get the cell types for the current cluster
    cell_type_names <- as.character(clusters_df_placenta_subset[clusters_df_placenta_subset$Cluster == i, ]$cell_type)
    
    # Count the occurrences of each cell type in the cluster
    cell_type_counts <- table(cell_type_names)
    
    # Get the top 3 most abundant cell types
    top3 <- names(sort(cell_type_counts, decreasing = TRUE)[1:3])
    # print(paste("Cluster: ", i))
    # print("Top 3 cell types:")
    # print(top3)
    
    # Check if the most abundant cell type is already used as a cluster ID
    chosen_cell_type <- top3[1]
    
    # Loop through the top 3 and select the first available one that hasn't been used yet
    for (cell_type in top3) {
        if (!(cell_type %in% used_cell_types)) {
            chosen_cell_type <- cell_type
            used_cell_types <- c(used_cell_types, chosen_cell_type)  # Mark this cell type as used
            break  # Exit the loop once a valid cell type is found
        }
    }
    
    # Append the chosen cell type to the new cluster ID list
    new.cluster.ids_placenta_subset_subset <- c(new.cluster.ids_placenta_subset_subset, chosen_cell_type)
}

new.cluster.ids_placenta_subset_subset




# Rename clusters according to most abundant cell type
new.cluster.ids_placenta_subset_subset <- c()
for (i in levels(clusters_df_placenta_subset$Cluster)){
    i <- as.numeric(i)
    cell_type_names <- as.character((clusters_df_placenta_subset[clusters_df_placenta_subset$Cluster == i,])$cell_type)
    cell_type_counts <- table(cell_type_names)
    print(paste("Cluster: ", i))
    most_abundant_cell_type <- names(as.list(head(sort(cell_type_counts, decreasing = TRUE), 1)))
    
    top3 <- names(as.list(head(sort(cell_type_counts, decreasing = TRUE), 3)))
    print(top3)
    print(most_abundant_cell_type)
    if (most_abundant_cell_type %in% clusters_df_placenta_subset){
        most_abundant_cell_type <- names(as.list(head(sort(cell_type_counts, decreasing = TRUE), 2)))[2]
        # most_abundant_cell_type <- gsub("^Placenta-", "", most_abundant_cell_type)
        new.cluster.ids_placenta_subset_subset <- c(new.cluster.ids_placenta_subset_subset, most_abundant_cell_type)
    } else{
        # most_abundant_cell_type <- gsub("^Placenta-", "", most_abundant_cell_type)
        new.cluster.ids_placenta_subset_subset <- c(new.cluster.ids_placenta_subset_subset, most_abundant_cell_type)
    }
    
    }

levels(clusters_df_placenta_subset$Cluster)

# Initialize the vector to store new cluster IDs
new.cluster.ids_placenta_subset_subset <- c()
used_cell_types <- c()  # Keep track of already assigned cell types

# Loop through each cluster
for (i in levels(clusters_df_placenta_subset$Cluster)){
    i <- as.numeric(i)
    
    # Get the cell types for the current cluster
    cell_type_names <- as.character(clusters_df_placenta_subset[clusters_df_placenta_subset$Cluster == i, ]$cell_type)
    
    # Count the occurrences of each cell type in the cluster
    cell_type_counts <- table(cell_type_names)
    
    # Get the top 3 most abundant cell types
    top3 <- names(sort(cell_type_counts, decreasing = TRUE)[1:3])
    # print(paste("Cluster: ", i))
    # print("Top 3 cell types:")
    # print(top3)
    
    # Check if the most abundant cell type is already used as a cluster ID
    chosen_cell_type <- top3[1]
    
    # Loop through the top 3 and select the first available one that hasn't been used yet
    for (cell_type in top3) {
        if (!(cell_type %in% used_cell_types)) {
            chosen_cell_type <- cell_type
            used_cell_types <- c(used_cell_types, chosen_cell_type)  # Mark this cell type as used
            break  # Exit the loop once a valid cell type is found
        }
    }
    
    # Append the chosen cell type to the new cluster ID list
    new.cluster.ids_placenta_subset_subset <- c(new.cluster.ids_placenta_subset_subset, chosen_cell_type)
}

new.cluster.ids_placenta_subset_subset


levels(placenta_subset)

names(new.cluster.ids_placenta_subset_subset) <- levels(placenta_subset)
placenta_subset <- RenameIdents(placenta_subset, new.cluster.ids_placenta_subset_subset)

umap_plot <- DimPlot(placenta_subset, reduction = "umap", label = TRUE, label.size = 4) + 
    NoLegend() +  ggtitle("UMAP Clustering of Placenta Data Subset") +
    labs(x = "UMAP 1", y = "UMAP 2")

umap_plot

ggsave("/mnt/DATA3/timo/plots/umap_plot_placenta_subset.pdf", plot = umap_plot, width = 10, height = 8)
