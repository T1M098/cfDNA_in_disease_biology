library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(uwot)
library(patchwork)
args <- commandArgs(TRUE)

set.seed(42)
data_raw <-  read.csv("/data/ranks/Ranks_FCC_RNA_cell_type_tissue.csv")
compartment_map = read.csv("/data/compartments/compartment_map.csv")

data <- data_raw

# Pivot the data to wide format
wide_data <- data %>%
    select(sample, cell_type_tissue, correlation) %>%  # Select relevant columns
    pivot_wider(names_from = sample, values_from = correlation, values_fill = list(correlation = 0))

head(wide_data)
wide_data <- wide_data %>% select(-cell_type_tissue) %>% as.matrix()
wide_data <- t(wide_data)
wide_data <- scale(wide_data)

# Perform UMAP analysis
umap_result <- umap(
    wide_data,
    n_neighbors = 15,  # Set based on your dataset
    min_dist = 0.1,    # Control the tightness of clusters
    metric = "correlation", # euclidean, manhattan
    n_components = 2,
    learning_rate = 1.0
)

# Convert UMAP results to a data frame
umap_df <- as.data.frame(umap_result)
colnames(umap_df) <- c("UMAP1", "UMAP2")

# Add the status column for coloring
# umap_df$status <- status


# add status or tumor fraction
metadata_df <- data %>% distinct(sample, status, tfx)
# head(metadata_df)

umap_df_plot <- umap_df
umap_df_plot$sample <- rownames(umap_df)
umap_df_plot <- left_join(umap_df_plot, metadata_df, by = "sample")
umap_df_plot <- umap_df_plot %>% select(order(names(umap_df_plot)))

# Visualize UMAP results with larger text elements
umap_status <- ggplot(umap_df_plot, aes(x = UMAP1, y = UMAP2, color = status)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "UMAP Clustering of Correlation Scores by Status",
    x = "UMAP1",
    y = "UMAP2",
    color = "Status"
  ) +
  theme(
    plot.title = element_text(size = 22),  # Increase title size
    axis.title = element_text(size = 16),  # Increase axis labels size
    axis.text = element_text(size = 14),   # Increase axis tick labels size
    legend.title = element_text(size = 16), # Increase legend title size
    legend.text = element_text(size = 14)   # Increase legend labels size
  )

umap_status

# ggsave("/plots/umap_plot_correlations_status.pdf", plot = umap_status, width = 10, height = 8)

# remove outliers tfx
umap_100 <- umap_df_plot %>% filter(0.0 <= tfx & tfx <= 1) # & status == "Healthy"

# Visualize UMAP results
umap_tfx_100 <- ggplot(umap_100, aes(x = UMAP1, y = UMAP2, color = tfx)) +
    geom_point(size = 3, alpha = 0.7) +
    theme_minimal() +
    labs(
        title = "Range 0-100%",
        x = "UMAP1",
        y = "UMAP2",
        color = "est. tfx"
    ) +
    scale_color_viridis_c(limits = c(min(umap_100 %>% distinct(tfx)), max(umap_100 %>% distinct(tfx)))) + 
    coord_cartesian(xlim = c(-4.5,4), ylim = c(-4.5,4)) +
  theme(
    plot.title = element_text(size = 22, hjust = 0.5),  # Increase title size
    axis.title = element_text(size = 16),  # Increase axis labels size
    axis.text = element_text(size = 14),   # Increase axis tick labels size
    legend.title = element_text(size = 16), # Increase legend title size
    legend.text = element_text(size = 14)   # Increase legend labels size
  )

umap_tfx_100

# remove outliers tfx
umap_2_50 <- umap_df_plot %>% filter(0.02 <= tfx & tfx <= .5) # & status == "Healthy"

# Visualize UMAP results
umap_tfx_2_50 <- ggplot(umap_2_50, aes(x = UMAP1, y = UMAP2, color = tfx)) +
    geom_point(size = 3, alpha = 0.7) +
    theme_minimal() +
    labs(
        title = "Range 2-50%",
        x = "UMAP1",
        y = "UMAP2",
        color = "est. tfx"
    ) +
    scale_color_viridis_c(limits = c(min(umap_2_50 %>% distinct(tfx)), max(umap_2_50 %>% distinct(tfx)))) + 
    coord_cartesian(xlim = c(-4.5,4), ylim = c(-4.5,4)) +
  theme(
    plot.title = element_text(size = 22, hjust = 0.5),  # Increase title size
    axis.title = element_text(size = 16),  # Increase axis labels size
    axis.text = element_text(size = 14),   # Increase axis tick labels size
    legend.title = element_text(size = 16), # Increase legend title size
    legend.text = element_text(size = 14)   # Increase legend labels size
  )

umap_tfx_2_50

# ggsave("/plots/umap_plot_correlations_tfx_2_to_50.pdf", plot = umap_tfx_2_50, width = 10, height = 8)

# Adjust spacing between the main title and subplot titles
combined_plot <- (umap_tfx_100 + umap_tfx_2_50) +
  plot_layout(nrow = 1, ncol = 2) +
  plot_annotation(
    title = "UMAP Clustering of Correlation Scores by est. Tumor Fraction",
    theme = theme(
      plot.title = element_text(size = 26, hjust = 0.5, margin = margin(b = 20)) # Increase bottom margin
    )
  )

# Display the plot
combined_plot

# ggsave("/plots/umap_plot_correlations_tfx_combined.pdf", plot = combined_plot, width = 14, height = 8)
