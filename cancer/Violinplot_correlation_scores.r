library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(uwot)
library(dplyr)
library(pheatmap)
library(viridis)
args <- commandArgs(TRUE)

set.seed(42)
data_rna_ctt <-  read.csv("data/ranks/Ranks_FCC_RNA_cell_type_tissue.csv")
data_rna_ct <-  read.csv("data/ranks/Ranks_FCC_RNA_cell_type.csv")
data_atac_ct <-  read.csv("data/ranks/Ranks_FCC_ATAC_cell_type.csv")
data_wps <- read.csv("data/ranks/Ranks_WPS_RNA_cell_type_tissue_correlations.tsv", sep="\t")
data_wps <- data_wps %>% rename(correlation = correlation_sample)

data_rna_ctt$type <- "FCC: RNA Cell Type Tissue"
data_rna_ct$type <- "FCC: RNA Cell Type"
data_atac_ct$type <- "FCC: ATAC-Seq Cell Type"
data_wps$type <- "WPS: RNA Cell Type Tissue"

data_rna_ctt <- data_rna_ctt %>% select(correlation, type)
data_rna_ct <- data_rna_ct %>% select(correlation, type)
data_atac_ct <- data_atac_ct %>% select(correlation, type)
data_wps <- data_wps %>% select(correlation, type)

merged_df <- rbind(data_rna_ctt, data_rna_ct, data_atac_ct, data_wps)
merged_df$type <- factor(merged_df$type, levels = c("FCC: ATAC-Seq Cell Type", "FCC: RNA Cell Type Tissue", "WPS: RNA Cell Type Tissue", "FCC: RNA Cell Type"))

# Calculate means and medians for each group
stats <- merged_df %>%
  group_by(type) %>%
  summarize(
    mean_correlation = mean(correlation),
    median_correlation = median(correlation)
  )

# Generate 4 discrete viridis colors
colors <- viridis(n = 4, direction = -1)

# Extract the two middle colors
middle_colors <- colors[2:3]

data_wps$type2 <- "WPS"
data_rna_ctt$type2 <- "FCC"

merged_fcc_wps <- rbind(data_wps, data_rna_ctt)
# Create the violin plot and add means and medians
ggplot(merged_fcc_wps, aes(x = type2, y = correlation, fill = type2)) +
  geom_violin(trim = FALSE, alpha = 0.7) +  # Violin plot
  geom_boxplot(width = 0.1, fill = "white") +  # Boxplot inside violin
  theme_minimal() +
  labs(title = "Distribution of PCC by Metric",
       x = "", 
       y = "PCC") +
  scale_fill_manual(values = c("#22A884", "#2A788E")) +
  theme(
    plot.title = element_text(size = 23, hjust = 0.5),  # Increase title font size
    axis.title.x = element_text(size = 20),  # Increase x-axis label font size
    axis.title.y = element_text(size = 20),  # Increase y-axis label font size
    legend.title = element_text(size = 20),  # Increase legend title font size
    legend.text = element_text(size = 17),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_blank()
  ) +
  guides(fill = guide_legend(title = "Metric"))  # Change legend title to "Metric"


# Save the plot to a file
# ggsave("/plots/violin_plot_WPS_vs_FCC.pdf", width = 8, height = 8)

# Viridis color palette
ggplot(merged_df, aes(x = type, y = correlation, fill = type)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white") +
  theme_minimal() +
  labs(title = "Distribution of PCC across Reference Data x Metric",
       x = "", 
       y = "PCC") +
  scale_fill_viridis(discrete = TRUE, direction = -1) +  # Use viridis palette
  theme(
    plot.title = element_text(size = 23, hjust = 0.5),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 17),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_blank()# element_text(angle = 45, hjust = 1)
  ) +
  guides(fill = guide_legend(title = "Combination"))

# Save the plot to a file
# ggsave("/plots/violin_plot_all_combinations.pdf", width = 12, height = 6)
