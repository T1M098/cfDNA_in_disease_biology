library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(pheatmap)
library(viridis)
custom_colors <- viridis(100, option = "viridis")[1:100]
library(patchwork)
args <- commandArgs(TRUE)

# Load data
atac <- read.csv("/data/evaluation_socres/FCC_ATAC_cell_type_eva_metrics_table.csv", header = TRUE)
rna_ct <- read.csv("//data/evaluation_socres/FCC_RNA_cell_type_eva_metrics_table.csv.csv", header = TRUE)
rna_ct_tissue <- read.csv("/data/evaluation_socres/FCC_RNA_cell_type_tissue_eva_metrics_table.csv", header = TRUE)
wps <- read.csv("/data/evaluation_socres/WPS_RNA_cell_type_tissue_eva_metrics_table.csv", header = TRUE)

# Combine all DataFrames
combined_df <- rbind(rna_ct_tissue, rna_ct, atac, wps)
# Reorder the levels of ref_data_type
combined_df$ref_data_type <- factor(combined_df$ref_data_type, levels = c("rna_ct_tissue", "rna", "atac", "WPS RNA Cell Type Tissue"))

combined_df$ref_data_type <- gsub("rna_ct_tissue", "FCC: RNA Cell Type Tissue", combined_df$ref_data_type)
combined_df$ref_data_type <- gsub("rna", "FCC: RNA Cell Type", combined_df$ref_data_type)
combined_df$ref_data_type <- gsub("atac", "FCC: ATAC-Seq Cell Type", combined_df$ref_data_type)
combined_df$ref_data_type <- gsub("WPS RNA Cell Type Tissue", "WPS: RNA Cell Type Tissue", combined_df$ref_data_type)

combined_df$CancerType <- gsub("Cancer", "C.", combined_df$CancerType)
combined_df$CancerType <- gsub("cancer", "C.", combined_df$CancerType)

# List of metrics
metrics <- c("AUC", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1_Score")

df_rna_ct_tissue <- combined_df %>% filter(ref_data_type == "FCC: RNA Cell Type Tissue")
heatmap_df_rna_ct_tissue <- df_rna_ct_tissue %>% select(-CancerType, -ref_data_type)

df <- df_rna_ct_tissue %>% select(-ref_data_type)

min_value_1 <- round(df %>%
  select(where(is.numeric)) %>%
  unlist() %>%
  min(na.rm = TRUE), 3)

# Reshape the data into long format
df_long <- df %>%
  pivot_longer(cols = -CancerType, names_to = "Metric", values_to = "Score")

# Define the desired custom order for CancerType
cancer_order <- c("Colorectal C.", "Ovarian C.", "Bile Duct C.", "Breast C.", "Gastric C.", "Lung C.", "Pancreatic C.")

# Reorder the CancerType factor in your dataset
df_long$CancerType <- factor(df_long$CancerType, levels = rev(cancer_order))

# Create the heatmap using ggplot2
heatmap_fcc_rna_ct_tissue <- ggplot(df_long, aes(x = Metric, y = CancerType, fill = Score)) +
  geom_tile() +
  scale_fill_viridis(option = "D", limits = c(min_value_1, 1)) +  # Color scale using viridis
  labs(title =  "Performance",
        subtitle = "FCC x RNA Cell Type Tissue",
        x = "Metric", y = "Observation") +
  theme_minimal() +
  theme(
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 18),  # Rotate x-axis labels for readability
        axis.text.y = element_text(size = 18),
        axis.title = element_blank(),  # Remove axis titles
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 22),
        plot.subtitle = element_text(size = 18, hjust = 0.5)
  )
heatmap_fcc_rna_ct_tissue

# Save the combined heatmaps as a PDF
# ggsave("/plots/eva_metrics_fcc_rna_ct_tissue_ordered.pdf", heatmap_fcc_rna_ct_tissue, width = 7, height = 7)

data_types <- c("FCC: RNA Cell Type Tissue", "FCC: RNA Cell Type", "FCC: ATAC-Seq Cell Type", "WPS: RNA Cell Type Tissue")

min_value <- round(combined_df %>%
  select(where(is.numeric)) %>%
  unlist() %>%
  min(na.rm = TRUE), 3)

# Initialize an empty list to store the heatmap plots
heatmaps <- list()

for (i in seq_along(data_types)) {
    data_type <- data_types[i]
    sub_df <- combined_df %>% filter(ref_data_type == data_type)
    df <- sub_df %>% select(-ref_data_type)

    # Reshape the data into long format
    df_long <- df %>%
        pivot_longer(cols = -CancerType, names_to = "Metric", values_to = "Score")

    # Define the desired custom order for CancerType
    cancer_order <- c("Colorectal C.", "Ovarian C.", "Bile Duct C.", "Breast C.", "Gastric C.", "Lung C.", "Pancreatic C.")
    
    # Reorder the CancerType factor in your dataset
    df_long$CancerType <- factor(df_long$CancerType, levels = rev(cancer_order))
    
    # Create the heatmap using ggplot2
    p <- ggplot(df_long, aes(x = Metric, y = CancerType, fill = Score)) +
        geom_tile() +
        scale_fill_viridis(option = "D", limits = c(min_value, 1)) +  # Color scale using viridis
        labs(title = gsub(":", "\n", data_type), x = "Metric", y = "Cancer Type") +
        theme_minimal() +
        theme(
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 16),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 18),  # Rotate x-axis labels for readability
            axis.text.y = element_text(size = 18),
            axis.title = element_blank(),  # Remove axis titles
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 22)# Remove grid lines
        )

    # Remove row labels for the first two heatmaps
    if (i == 2 || i == 3 || i == 4) {
        p <- p + theme(axis.text.y = element_blank())  # Remove y-axis text (row labels)
    }

    # Add the legend only for the third heatmap
    if (i == 4) {
        p <- p + guides(fill = guide_colorbar(title = "Score"))  # Show legend for the third heatmap
    } else {
        p <- p + theme(legend.position = "none")  # Remove legend for the first two heatmaps
    }

    # Add the ggplot object to the list
    heatmaps[[data_type]] <- p
}

# Combine all heatmaps into a single plot layout with wrap_plots
combined_heatmaps <- wrap_plots(heatmaps, ncol = 4)
combined_heatmaps


# Save the combined heatmaps as a PDF
# ggsave("/plots/eva_metrics_combined_heatmaps_per_ref_data_wps_ordered.pdf", combined_heatmaps, width = 16, height = 6)

# Load data
pcc_df <- read.csv("/main/body/outputs/top_highest_lowest_pcc_values.csv", header = TRUE)
svm_df <- read.csv("/main/body/outputs/feature_weights_loosvm_all_cancer_types.csv", header = TRUE)
svm_df <- svm_df %>% rename(cell_type = Feature)

# Rank by PCC
pcc_df <- pcc_df %>%
    group_by(CancerType) %>%
    arrange(desc(PCC)) %>%  # Use `asc(PCC)` for ascending order if needed
    mutate(rank_score_PCC = row_number()) %>%
    arrange(CancerType)
    
# Rank by feature weight (AvgWeight)
svm_df <- svm_df %>%
    group_by(CancerType) %>%
    arrange(desc(AvgWeight)) %>%  # Use `asc(AvgWeight)` for ascending order if needed
    mutate(rank_AvgWeight = row_number()) %>%
    arrange(CancerType)

n <- 50
# Select the top n highest and lowest for PCC
top_n_pcc_highest <- pcc_df %>%
    group_by(CancerType) %>%
    top_n(n, PCC) %>%
    filter(p_value < 0.05)

top_n_pcc_lowest <- pcc_df %>%
    group_by(CancerType) %>%
    top_n(-n, PCC) %>%
    filter(p_value < 0.05)

# Select the top n highest and lowest for Feature Weights
top_n_svm_highest <- svm_df %>%
    group_by(CancerType) %>%
    top_n(n, AvgWeight)

top_n_svm_lowest <- svm_df %>%
    group_by(CancerType) %>%
    top_n(-n, AvgWeight)

# Transform ranks for top n lowest, 447 rank becomes 1
top_n_pcc_lowest <- top_n_pcc_lowest %>%
  mutate(rank_score_PCC_new = 448 - rank_score_PCC)

top_n_svm_lowest <- top_n_svm_lowest %>%
  mutate(rank_AvgWeight_new = 448 - rank_AvgWeight)

# Compare PCC and feature weight top 25 highest
common_highest <- inner_join(top_n_pcc_highest, top_n_svm_highest, by = c("CancerType", "cell_type"))

# Compare PCC and feature weight top 25 lowest
common_lowest <- inner_join(top_n_pcc_lowest, top_n_svm_lowest, by = c("CancerType", "cell_type"))

common_highest_save <- common_highest %>% select(CancerType, cell_type, rank_AvgWeight, rank_score_PCC, AvgWeight, PCC, p_value)
common_lowest_save <- common_lowest %>% select(CancerType, cell_type, rank_AvgWeight, rank_score_PCC, AvgWeight, PCC, p_value)

common_highest_save$cell_type <- gsub("_", " ", common_highest_save$cell_type)
common_lowest_save$cell_type <- gsub("_", " ", common_lowest_save$cell_type)

common_highest_save <- common_highest_save %>%
  mutate_at(vars(AvgWeight, PCC), round, 3)

common_lowest_save <- common_lowest_save %>%
  mutate_at(vars(AvgWeight, PCC), round, 3)

common_highest_save$p_value <- format(common_highest_save$p_value, scientific = TRUE, digits = 3)
common_lowest_save$p_value <- format(common_lowest_save$p_value, scientific = TRUE, digits = 3)

colnames(common_highest_save) <- c("Cancer Type", "Cell Type", "Rank AvgWeight", "Rank PCC Score", 
                   "Avg Weight", "PCC Score", "p-value PCC Score")
colnames(common_lowest_save) <- c("Cancer Type", "Cell Type", "Rank AvgWeight", "Rank PCC Score", 
                   "Avg Weight", "PCC Score", "p-value PCC Score")

common_highest_save[] <- lapply(common_highest_save, function(x) if(is.character(x)) gsub('^"|"$', '', x) else x)
common_lowest_save[] <- lapply(common_lowest_save, function(x) if(is.character(x)) gsub('^"|"$', '', x) else x)

common_lowest_save

# write.csv(common_highest_save, "/outputs/highest_feature_weights_pcc_scores_50.csv", row.names = FALSE)
# write.csv(common_lowest_save, "/outputs/lowest_feature_weights_pcc_scores_50.csv", row.names = FALSE)
