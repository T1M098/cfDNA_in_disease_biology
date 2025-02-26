library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(pheatmap)
library(viridis)
custom_colors <- viridis(100, option = "viridis")[1:100]
library(patchwork)

# Load data
atac <- read.csv("/data/evaluation_scores/FCC_ATAC_cell_type_eva_metrics_table.csv", header = TRUE)
rna_ct <- read.csv("/data/evaluation_scores/FCC_RNA_cell_type_eva_metrics_table.csv", header = TRUE)
rna_ct_tissue <- read.csv("/data/evaluation_scores/FCC_RNA_cell_type_tissue_eva_metrics_table.csv", header = TRUE)
ovr <- read.csv("/data/evaluation_scores/FCC_RNA_cell_type_tissue_OvR_eva_metrics_table.csv.csv", header = TRUE)
wps <- read.csv("/data/evaluation_scores/WPS_RNA_cell_type_tissue_eva_metrics_table.csv", header = TRUE)

# Combine all DataFrames
combined_df <- rbind(rna_ct_tissue, rna_ct, atac, wps, ovr)
# Reorder the levels of ref_data_type
combined_df$ref_data_type <- factor(combined_df$ref_data_type, levels = c("rna_ct_tissue", "rna", "atac", "WPS RNA Cell Type Tissue", "rna_ct_tissue_ovr"))

combined_df$ref_data_type <- gsub("rna_ct_tissue", "FCC: RNA Cell Type Tissue", combined_df$ref_data_type)
combined_df$ref_data_type <- gsub("rna", "FCC: RNA Cell Type", combined_df$ref_data_type)
combined_df$ref_data_type <- gsub("atac", "FCC: ATAC-Seq Cell Type", combined_df$ref_data_type)
combined_df$ref_data_type <- gsub("WPS RNA Cell Type Tissue", "WPS: RNA Cell Type Tissue", combined_df$ref_data_type)

combined_df$CancerType <- gsub("Cancer", "C.", combined_df$CancerType)
combined_df$CancerType <- gsub("cancer", "C.", combined_df$CancerType)

# List of metrics
metrics <- c("AUC", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1_Score")

head(combined_df)

ovr <- combined_df %>% filter(ref_data_type == "FCC: RNA Cell Type Tissue_ovr")

df_rna_ct_tissue <- combined_df %>% filter(ref_data_type == "FCC: RNA Cell Type Tissue")
heatmap_df_rna_ct_tissue <- df_rna_ct_tissue %>% select(-CancerType, -ref_data_type)

# Assuming your df_rna_ct_tissue dataframe is already defined
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

# Create a new column for text color based on the score value
df_long$text_color <- ifelse(df_long$Score < 0.5, "grey", "black")

# Create the heatmap using ggplot2 with scores inside the tiles
heatmap_fcc_rna_ct_tissue_score <- ggplot(df_long, aes(x = Metric, y = CancerType, fill = Score)) +
  geom_tile() +
  geom_text(aes(label = round(Score, 2), color = text_color), size = 4) +  # Adding text inside heatmap cells
  scale_fill_viridis(option = "D", limits = c(min_value_1, 1)) +  # Color scale using viridis
  scale_color_identity() +  # Apply color without changing scale
  labs(title =  "Performance SVM",
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

# Display the heatmap
heatmap_fcc_rna_ct_tissue_score


# Save the combined heatmaps as a PDF
# ggsave("/plots/eva_metrics_fcc_rna_ct_tissue_ordered_scores.pdf", heatmap_fcc_rna_ct_tissue_score, width = 7, height = 7)

# Assuming your df_rna_ct_tissue dataframe is already defined
df_ovr <- ovr %>% select(-ref_data_type)

min_value_1 <- round(df_ovr %>%
  select(where(is.numeric)) %>%
  unlist() %>%
  min(na.rm = TRUE), 3)

# Reshape the data into long format
df_long <- df_ovr %>%
  pivot_longer(cols = -CancerType, names_to = "Metric", values_to = "Score")

# Define the desired custom order for CancerType
cancer_order <- c("Colorectal C.", "Ovarian C.", "Bile Duct C.", "Breast C.", "Gastric C.", "Lung C.", "Pancreatic C.")

# Reorder the CancerType factor in your dataset
df_long$CancerType <- factor(df_long$CancerType, levels = rev(cancer_order))

# Create a new column for text color based on the score value
df_long$text_color <- ifelse(df_long$Score < 0.5, "grey", "black")

# Create the heatmap using ggplot2 with scores inside the tiles
heatmap_fcc_ovr_score <- ggplot(df_long, aes(x = Metric, y = CancerType, fill = Score)) +
  geom_tile() +
  geom_text(aes(label = round(Score, 2), color = text_color), size = 4) +  # Adding text inside heatmap cells
  scale_fill_viridis(option = "D", limits = c(min_value_1, 1)) +  # Color scale using viridis
  scale_color_identity() +  # Apply color without changing scale
  labs(title =  "Performance OvR SVM",
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

# Display the heatmap
heatmap_fcc_ovr_score

# Save the combined heatmaps as a PDF
# ggsave("/plots/eva_metrics_fcc_rna_ct_tissue_ordered_scores_ovr.pdf", heatmap_fcc_ovr_score, width = 7, height = 7)

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

    # Create a new column for text color based on the score value
    df_long$text_color <- ifelse(df_long$Score < 0.3, "grey", "black")

    # Create the heatmap using ggplot2
    p <- ggplot(df_long, aes(x = Metric, y = CancerType, fill = Score)) +
        geom_tile() +
        geom_text(aes(label = round(Score, 2), color = text_color), size = 4) +  # Adding text inside heatmap cells
        scale_fill_viridis(option = "D", limits = c(min_value, 1)) +  # Color scale using viridis
        scale_color_identity() +  # Apply color without changing scale
        labs(title = gsub(":", "\n", data_type), x = "Metric", y = "Cancer Type") +
        theme_minimal() +
        theme(
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 16),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 18),  # Rotate x-axis labels for readability
            axis.text.y = element_text(size = 18),
            axis.title = element_blank(),  # Remove axis titles
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 22)  # Remove grid lines
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
combined_heatmaps_scores <- wrap_plots(heatmaps, ncol = 4)
combined_heatmaps_scores


# Save the combined heatmaps as a PDF
# ggsave("/plots/eva_metrics_combined_heatmaps_per_ref_data_wps_ordered_scores.pdf", combined_heatmaps_scores, width = 16, height = 6)
