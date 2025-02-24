library(data.table)
library(Rtsne)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(igraph)
library(FNN)
library(cowplot)
library(tidyverse)

library(viridis)
# Extract a color from the 'viridis' palette
viridis(1, option = "viridis")


data <- read.csv("/data/ranks/Ranks_FCC_RNA_cell_type_tissue.csv", header = TRUE)
compartment_map <- read.csv("/data/compartments/compartment_map.csv")
# compartment_map <- compartment_map %>% mutate(cell_type_tissue = str_replace(cell_type_tissue, ",", ""))
# compartment_map <- compartment_map %>% mutate(cell_type_tissue = str_replace_all(cell_type_tissue, "-", " "))
# compartment_map$cell_type_tissue <- str_trim(compartment_map$cell_type_tissue)
# write.csv(compartment_map, file = "/mnt/DATA3/timo/data/compartment_map.csv", row.names = FALSE)

combined <- data %>% left_join(compartment_map, by = c("cell_type_tissue" = "cell_type_tissue"))

healthy <- combined %>% filter(status == "Healthy")
healthy <- healthy %>% mutate(tissue_type = str_extract(cell_type_tissue, "(?<=_).*$"))

# Compute the median rank for each tissue type
healthy_medians_compartment <- healthy %>%
  group_by(compartment) %>%
  summarize(median_rank = median(rank, na.rm = TRUE), .groups = "drop")

# Merge median ranks back into the original data
healthy_compartment <- healthy %>%
  left_join(healthy_medians_compartment, by = "compartment")

# Create the boxplot with colors mapped to the median rank
boxplot_rank_compartment <- ggplot(healthy_compartment, aes(x = reorder(compartment, -rank, FUN = median), y = rank, fill = median_rank)) +
  geom_boxplot(color = "black") +  
  theme_minimal() +
  labs(
    title = "Rank Distribution across Compartments of Healthy Samples",
    x = "",
    y = "Rank",
    fill = "Median Rank"
  ) +
  scale_y_continuous(trans = 'reverse') +  # Reverse y-axis to match your example
  coord_flip() +  # Flip coordinates to make it horizontal
  scale_fill_viridis(option = "D", direction = -1) +  # Apply Viridis colormap with limits
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),  # Center and increase title size
    axis.text.x = element_text(size = 18),  # Increase x-axis label size
    axis.text.y = element_text(size = 18),  # Increase y-axis label size
    axis.title.x = element_text(size = 20),
    legend.position = "none"  # Keep legend on the right
  )

# Display the plot
print(boxplot_rank_compartment)


# ggsave("/plots/rank_distribution_compartment.pdf", plot = boxplot_rank_compartment, width = 9, height = 5)

# Compute the median rank for each tissue type
healthy_medians_tissue <- healthy %>%
  group_by(tissue_type) %>%
  summarize(median_rank = median(rank, na.rm = TRUE), .groups = "drop")

# Merge median ranks back into the original data
healthy_tissue <- healthy %>%
  left_join(healthy_medians_tissue, by = "tissue_type")

# Create the boxplot with colors mapped to the median rank
boxplot_tissue <- ggplot(healthy_tissue, aes(x = reorder(tissue_type, -rank, FUN = median), y = rank, fill = median_rank)) +
  geom_boxplot(color = "black") +  
  theme_minimal() +
  labs(
    title = "Rank Distribution across Tissue Types of Healthy Samples",
    x = "",
    y = "Rank",
    fill = "Median Rank"
  ) +
  scale_y_continuous(trans = 'reverse') +  # Reverse y-axis to match your example
  coord_flip() +  # Flip coordinates to make it horizontal
  scale_fill_viridis(option = "D", direction = -1) +  # Apply Viridis colormap with limits
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),  # Center and increase title size
    axis.text.x = element_text(size = 18),  # Increase x-axis label size
    axis.text.y = element_text(size = 18),  # Increase y-axis label size
    axis.title.x = element_text(size = 20),  # Increase y-axis title size
    legend.position = "none"  # Keep legend on the right
  )

# Display the plot
print(boxplot_tissue)


# ggsave("/plots/rank_distribution_tissue_type.pdf", plot = boxplot_tissue, width = 10, height = 10)

healthy_plasma_cell <- healthy[grepl("plasma cell", healthy$cell_type, ignore.case = TRUE), ]

# Compute the median rank for each tissue type
healthy_medians_plasma_cell <- healthy_plasma_cell %>% 
  group_by(tissue_type) %>%
  summarize(median_rank = median(rank, na.rm = TRUE), .groups = "drop")

# Merge median ranks back into the original data
healthy_tissue_plasma_cell <- healthy_plasma_cell %>%
  left_join(healthy_medians_plasma_cell, by = "tissue_type")

# Create the boxplot with colors mapped to the median rank
boxplot_tissue_plasma_cells <- ggplot(healthy_tissue_plasma_cell, aes(x = reorder(tissue_type, -rank, FUN = median), y = rank, fill = median_rank)) +
  geom_boxplot(color = "black") +  
  theme_minimal() +
  labs(
    title = "Rank Distribution of Plasma Cells across Tissue Types of Healthy Samples",
    x = "",
    y = "Rank",
    fill = "Median Rank"
  ) +
  scale_y_continuous(trans = 'reverse') +  # Reverse y-axis to match your example
  coord_flip() +  # Flip coordinates to make it horizontal
  scale_fill_viridis(option = "D", direction = -1) +  # Apply Viridis colormap with limits
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),  # Center and increase title size
    axis.text.x = element_text(size = 18),  # Increase x-axis label size
    axis.text.y = element_text(size = 18),  # Increase y-axis label size
    axis.title.x = element_text(size = 20),  # Increase y-axis title size
    legend.position = "none"  # Keep legend on the right
  )

# Display the plot
print(boxplot_tissue_plasma_cells)

# ggsave("/plots/rank_distribution_plasma_cells_tissue_type.pdf", plot = boxplot_tissue_plasma_cells, width = 12, height = 12)
