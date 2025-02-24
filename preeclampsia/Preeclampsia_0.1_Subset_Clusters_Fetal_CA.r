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

# install.packages("devtools")
# devtools::install_github("immunogenomics/presto")

# packageVersion("presto")

# packageVersion("Seurat")

print("Loading .rds test file.")
fetal_cell_atlas <- readRDS("/mnt/DATA3/timo/data/placenta/fetal_cell_atlas_placenta_seurat_object_norm.rds")
print("File loaded.")

(fetal_cell_atlas@assays)

unique(fetal_cell_atlas@meta.data$cell_type)

fetal_cell_atlas <- FindVariableFeatures(fetal_cell_atlas, selection.method = "vst", nfeatures = 2000)

DefaultAssay(fetal_cell_atlas) <- "RNA"

fetal_cell_atlas <- ScaleData(fetal_cell_atlas, features = rownames(fetal_cell_atlas))

fetal_cell_atlas <- RunPCA(fetal_cell_atlas, features = VariableFeatures(fetal_cell_atlas), dims = 1:20)

# Examine and visualize PCA results a few different ways
print(fetal_cell_atlas[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(fetal_cell_atlas, reduction = "pca") + NoLegend()

fetal_cell_atlas <- FindNeighbors(fetal_cell_atlas, dims = 1:10)
fetal_cell_atlas <- FindClusters(fetal_cell_atlas, resolution = 0.1)

fetal_cell_atlas <- RunUMAP(fetal_cell_atlas, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(fetal_cell_atlas, reduction = "umap", label = TRUE) + ggtitle("UMAP of original Fetal CA placenta cells Clusters")

saveRDS(fetal_cell_atlas, file = "/mnt/DATA3/timo/data/placenta/fetal_cell_atlas_original_clusters.rds")

print("Loading .rds test file.")
fetal_cell_atlas <- readRDS("/mnt/DATA3/timo/data/placenta/fetal_cell_atlas_original_clusters.rds")
print("File loaded.")

DimPlot(fetal_cell_atlas, reduction = "umap", label = TRUE) + ggtitle("UMAP of original Fetal CA placenta cells Clusters")

# merge cluster 7 and 10 (most abundant cell type = Placenta-Syncytiotrophoblasts and villous cytotrophoblasts) with cluster 0 (most abundant cell type = Placenta-Syncytiotrophoblasts and villous cytotrophoblasts)
Idents(fetal_cell_atlas)[Idents(fetal_cell_atlas) %in% c("7", "10")] <- "0"

fetal_cell_atlas@meta.data$seurat_clusters_new <- fetal_cell_atlas@meta.data$seurat_clusters
fetal_cell_atlas@meta.data$seurat_clusters_new <- as.character(fetal_cell_atlas@meta.data$seurat_clusters_new)  # Ensure it's character for modification
fetal_cell_atlas@meta.data$seurat_clusters_new[fetal_cell_atlas@meta.data$seurat_clusters_new %in% c("7", "10")] <- "0"
fetal_cell_atlas@meta.data$seurat_clusters_new <- as.factor(fetal_cell_atlas@meta.data$seurat_clusters_new)  # Convert back to factor if necessary

clusters_df_fetal_cell_atlas <- data.frame(Cell = rownames(fetal_cell_atlas@meta.data),
                          cell_type = fetal_cell_atlas@meta.data$cell_type,
                          Cluster = fetal_cell_atlas$seurat_clusters_new)

levels(clusters_df_fetal_cell_atlas$Cluster)

new.cluster.ids_fetal_cell_atlas <- c()
for (i in levels(clusters_df_fetal_cell_atlas$Cluster)){
    i <- as.numeric(i)
    cell_type_names <- as.character((clusters_df_fetal_cell_atlas[clusters_df_fetal_cell_atlas$Cluster == i,])$cell_type)
    cell_type_counts <- table(cell_type_names)
    print(paste("Cluster: ", i))
    most_abundant_cell_type <- names(as.list(head(sort(cell_type_counts, decreasing = TRUE), 1)))
    
    top3 <- names(as.list(head(sort(cell_type_counts, decreasing = TRUE), 3)))
    print(most_abundant_cell_type)
    if (most_abundant_cell_type %in% clusters_df_fetal_cell_atlas){
        most_abundant_cell_type <- names(as.list(head(sort(cell_type_counts, decreasing = TRUE), 2)))[2]
        most_abundant_cell_type <- gsub("^Placenta-", "", most_abundant_cell_type)
        new.cluster.ids_fetal_cell_atlas <- c(new.cluster.ids_fetal_cell_atlas, most_abundant_cell_type)
    } else{
        most_abundant_cell_type <- gsub("^Placenta-", "", most_abundant_cell_type)
        new.cluster.ids_fetal_cell_atlas <- c(new.cluster.ids_fetal_cell_atlas, most_abundant_cell_type)
    }
    
    }

new.cluster.ids_fetal_cell_atlas

for (i in 0:(length(unique(clusters_df_fetal_cell_atlas$Cluster))-1)){
    cell_types <- as.character((clusters_df_fetal_cell_atlas[clusters_df_fetal_cell_atlas$Cluster == i,])$cell_type)
    sorted_cell_types <- sort(table(cell_types), decreasing = TRUE)
    print(i)
    print(sorted_cell_types)
    }

(new.cluster.ids_fetal_cell_atlas)

levels(fetal_cell_atlas)

(new.cluster.ids_fetal_cell_atlas)

names(new.cluster.ids_fetal_cell_atlas) <- levels(fetal_cell_atlas)
fetal_cell_atlas <- RenameIdents(fetal_cell_atlas, new.cluster.ids_fetal_cell_atlas)

DimPlot(fetal_cell_atlas, reduction = "umap", label = TRUE, label.size = 3) + NoLegend() + ggtitle("UMAP of merged Fetal CA Clusters")

#### sample 2k cells

s <- 0
for (cluster_id in unique(Idents(fetal_cell_atlas))) {

    cells_in_cluster <- WhichCells(fetal_cell_atlas, idents = cluster_id)

    n_cells <- as.integer(length(cells_in_cluster))
    print(n_cells)
    print(cluster_id)
    s <- s + n_cells
}
s

sampled_cells <- list()
sample_size <- 3000
cluster_ids <- unique(Idents(fetal_cell_atlas))
for (cluster_id in cluster_ids) {
    cells_in_cluster <- WhichCells(fetal_cell_atlas, idents = cluster_id)
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
fetal_cell_atlas_subset <- subset(fetal_cell_atlas, cells = cells_to_keep)

fetal_cell_atlas_subset <- FindVariableFeatures(fetal_cell_atlas_subset, selection.method = "vst", nfeatures = 2000)

fetal_cell_atlas_subset <- ScaleData(fetal_cell_atlas_subset, features = rownames(fetal_cell_atlas_subset))
fetal_cell_atlas_subset <- RunPCA(fetal_cell_atlas_subset, features = VariableFeatures(fetal_cell_atlas_subset), dims = 1:20)

DimPlot(fetal_cell_atlas_subset, reduction = "pca") + NoLegend()

fetal_cell_atlas_subset <- FindNeighbors(fetal_cell_atlas_subset, dims = 1:10)
fetal_cell_atlas_subset <- FindClusters(fetal_cell_atlas_subset, resolution = 0.1)

fetal_cell_atlas_subset <- RunUMAP(fetal_cell_atlas_subset, dims = 1:10)

clusters_df_fetal_cell_atlas_subset <- data.frame(Cell = rownames(fetal_cell_atlas_subset@meta.data),
                          cell_type = fetal_cell_atlas_subset@meta.data$cell_type,
                          Cluster = fetal_cell_atlas_subset$seurat_clusters_new)

# Rename clusters according to most abundant cell type
new.cluster.ids_fetal_cell_atlas_subset <- c()
for (i in levels(clusters_df_fetal_cell_atlas_subset$Cluster)){
    i <- as.numeric(i)
    cell_type_names <- as.character((clusters_df_fetal_cell_atlas_subset[clusters_df_fetal_cell_atlas_subset$Cluster == i,])$cell_type)
    cell_type_counts <- table(cell_type_names)
    print(paste("Cluster: ", i))
    most_abundant_cell_type <- names(as.list(head(sort(cell_type_counts, decreasing = TRUE), 1)))
    
    top3 <- names(as.list(head(sort(cell_type_counts, decreasing = TRUE), 3)))
    print(most_abundant_cell_type)
    if (most_abundant_cell_type %in% clusters_df_fetal_cell_atlas){
        most_abundant_cell_type <- names(as.list(head(sort(cell_type_counts, decreasing = TRUE), 2)))[2]
        most_abundant_cell_type <- gsub("^Placenta-", "", most_abundant_cell_type)
        new.cluster.ids_fetal_cell_atlas_subset <- c(new.cluster.ids_fetal_cell_atlas_subset, most_abundant_cell_type)
    } else{
        most_abundant_cell_type <- gsub("^Placenta-", "", most_abundant_cell_type)
        new.cluster.ids_fetal_cell_atlas_subset <- c(new.cluster.ids_fetal_cell_atlas_subset, most_abundant_cell_type)
    }
    
    }

new.cluster.ids_fetal_cell_atlas_subset

levels(fetal_cell_atlas_subset)

names(new.cluster.ids_fetal_cell_atlas_subset) <- levels(fetal_cell_atlas_subset)
fetal_cell_atlas_subset <- RenameIdents(fetal_cell_atlas_subset, new.cluster.ids_fetal_cell_atlas_subset)

DimPlot(fetal_cell_atlas_subset, reduction = "umap", label = TRUE, label.size = 3) + NoLegend() + ggtitle("UMAP of subset of Fetal CA merged Cluster")

saveRDS(fetal_cell_atlas_subset, file = "/mnt/DATA3/timo/data/placenta/fetal_CA_subset3k_merged_clusters.rds")

fetal_cell_atlas_subset <- readRDS("/mnt/DATA3/timo/data/placenta/fetal_CA_subset3k_merged_clusters.rds")

umap_plot <- DimPlot(fetal_cell_atlas_subset, reduction = "umap", label = TRUE, label.size = 4) + 
    NoLegend() + ggtitle("UMAP Clustering of Fetal Cell Atlas Subset") +
    scale_x_continuous(expand = expansion(add = c(3, 3))) +
    labs(x = "UMAP 1", y = "UMAP 2")
umap_plot

ggsave("/mnt/DATA3/timo/plots/umap_plot_fetal_cell_atlas_subset.pdf", plot = umap_plot, width = 10, height = 8)
