rm(list = ls())

# Used to understand and modify functions.R which includes wilcox tests for volcano plots and walktrap algorithm
library(data.table)
library(Rtsne)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(igraph)
library(FNN)
library(cowplot)
library(tidyverse)
library(rstatix)
library(gridExtra)
library(viridis)
library(patchwork)
functions <- file.path("//stanley/case_control/functions.R")
source(functions)

data <- read.csv("data/ranks/Ranks_FCC_RNA_cell_type_tissue.csv", header = TRUE)
# data <- data %>% rename(cell_type_tissue = cell_type)
# data$cell_type_tissue <- gsub("^RNA\\.", "", data$cell_type_tissue)
# # Replace all '.' or '..' with spaces
# data$cell_type_tissue <- gsub("_", " ", data$cell_type_tissue)  # Replace hyphens with spaces
# data$cell_type_tissue <- gsub("\\.{1,2}", " ", data$cell_type_tissue)
# data$cell_type_tissue <- gsub("\\s{2,}", " ", data$cell_type_tissue) # replace double spaces

# # # Replace the last space with an underscore in the 'cell_type_tissue' column
# data$cell_type_tissue <- sub(" (?!.* )", "_", data$cell_type_tissue, perl = TRUE)
# data$cell_type_tissue <- gsub(" Small_", "_Small ", data$cell_type_tissue)
# data$cell_type_tissue <- gsub(" Large_", "_Large ", data$cell_type_tissue)
# data$cell_type_tissue <- gsub(" Salivary_", "_Salivary ", data$cell_type_tissue)
# data$cell_type_tissue <- gsub(" Bone_", "_Bone ", data$cell_type_tissue)
# data$cell_type_tissue <- gsub(" Lymph_", "_Lymph ", data$cell_type_tissue)
# # Save your dataframe as a CSV file
# data <- data %>%
#   filter(!(grepl("Prostate", cell_type_tissue) & gender == "F"))
# write.csv(data, file = "/mnt/DATA3/timo/centers/main/body/outputs/combined_518_FCC_final.csv", row.names = FALSE)
# data
compartment_map <- read.csv("/data/compartments/compartment_map.csv")
# remove "," and "-" from cell_type for immune cells, e.g. effector CD8-positive, alpha-beta T cell
compartment_map <- compartment_map %>% mutate(cell_type_tissue = str_replace(cell_type_tissue, ",", ""))
compartment_map <- compartment_map %>% mutate(cell_type_tissue = str_replace_all(cell_type_tissue, "-", " "))
compartment_map$cell_type_tissue <- str_trim(compartment_map$cell_type_tissue)
# compartment_map %>% distinct(cell_type_tissue)

test <- data %>% left_join(compartment_map, by = c("cell_type_tissue" = "cell_type_tissue"))
# test %>% filter(if_any(everything(), is.na)) %>% distinct(cell_type)

crc_test <- data %>% filter(status %in% c("Colorectal Cancer", "Healthy"))

plotDifferentialCellTypes <- function(data, cancer_type, output_dir) {
    
    allcomp <- c()
    cell <- unique(data$cell_type_tissue)
    
    #run two-sided wilcox test and compute mean foldchange for each cell type
    for(i in 1:length(cell)) {
    df <- subset(data, cell_type_tissue %in% cell[i])
    pheno.counts <- as.data.table(table(df$status))
    labels <- paste0(pheno.counts[1,1],".n.", pheno.counts[1,2], ".vrs.", pheno.counts[2,1], ".n.", pheno.counts[2,2])
    df.wilcox <- wilcox_test(df, rank ~ status, ref.group = pheno.counts[2,1]$V1)
    df.mean <- df %>% group_by(status) %>% summarize(Mean = mean(rank, na.rm=TRUE)) %>% as.data.frame()
    df.mean.FC <- df.mean$Mean[2]/df.mean$Mean[1] %>% as.data.frame()
    colnames(df.mean.FC) <- c("foldchange")
    df.mean.FC$cell_type_tissue <- cell[i]
    df.all <- cbind(df.wilcox, df.mean.FC)
    allcomp <- rbind(allcomp, df.all)
    }
    
    #multiple testing correction
    allcomp$p.adj <- p.adjust(allcomp$p, method="fdr")

    #order by p.adj
    allcomp <- allcomp[order(allcomp[,10]),]
    
    #create color for volcano plot
    allcomp$diffranked <- "Not Sig"
    allcomp$diffranked[allcomp$foldchange > 1 & allcomp$p.adj < 0.05] <- "Up"
    allcomp$diffranked[allcomp$foldchange < 1 & allcomp$p.adj < 0.05] <- "Down"

    allcomp$delabel <- NA
    allcomp[allcomp$diffranked == "Down" | allcomp$diffranked == "Up",]$delabel <- allcomp[allcomp$diffranked == "Down" | allcomp$diffranked == "Up",]$cell_type

    # add compartment and cell_type_tissue
    allcomp$cell_type_tissue <- str_trim(allcomp$cell_type_tissue)
    allcomp <- allcomp %>% left_join(compartment_map, by = c("cell_type_tissue" = "cell_type_tissue"))
    allcomp <- allcomp %>%
    # Create a new column 'top_label' based on the conditions
    mutate(
    top_label = case_when(
      diffranked == "Down" & foldchange %in% head(sort(foldchange[diffranked == "Down"], decreasing = FALSE), 10) ~ cell_type_tissue,  # Top 10 lowest foldchange for "Down"
      diffranked == "Up" & foldchange %in% head(sort(foldchange[diffranked == "Up"], decreasing = TRUE), 10) ~ cell_type_tissue,    # Top 10 highest foldchange for "Up"
      TRUE ~ NA_character_  # Set NA for all other rows
    )
    )

    line_up_fc_value <- min(allcomp %>% filter(diffranked =="Up", !is.na(top_label)) %>% pull(foldchange))
    line_down_fc_value <- max(allcomp%>% filter(diffranked =="Down", !is.na(top_label)) %>% pull(foldchange))
    
    p_sgnf_base_line <- min(allcomp %>% filter(diffranked == "Not Sig") %>% pull(p))

    xlim_up <- max(allcomp %>% pull(foldchange))
    xlim_up <- xlim_up * 1.3

    xlim_down <- min(allcomp %>% pull(foldchange))
    xlim_down <- xlim_down * 0.7
    
    #set colors for volcano plot
    mycolors <- c("#541352FF", "#10a53dFF", "grey80")
    names(mycolors) <- c("Down", "Up", "Not Sig")

    my_compartment_shapes <- c(
      "immune" = 15,       
      "endothelial" = 16,  
      "epithelial" = 17,  
      "stromal" = 3,      
      "unknown" = 18       
    )

    plot <- ggplot(data=allcomp, aes(x=log2(foldchange), y=-log10(p), color=diffranked, shape=compartment, label=top_label)) +
        geom_vline(xintercept = log2(line_up_fc_value), color = "gray", linetype = "dotted", size = 0.2) +   
        geom_vline(xintercept = log2(line_down_fc_value), color = "gray", linetype = "dotted", size = 0.2) +
        geom_hline(yintercept = -log10(p_sgnf_base_line), color = "gray", linetype = "dotted", size = 0.2) +
        
        geom_point() +
        geom_label_repel(aes(label = ifelse(log2(foldchange) > log2(line_down_fc_value) | 
                                          log2(foldchange) < log2(line_up_fc_value), top_label, "")), size = 3, show.legend = FALSE) +
        
        scale_color_manual(values = mycolors) +
        scale_shape_manual(values = my_compartment_shapes) +
    
        guides(
            shape = guide_legend(order = 2, size = 3),
            color = guide_legend(order = 1, override.aes = list(shape = 16, size = 3))  
        ) +
    
        theme_classic() +
        ggtitle(paste(cancer_type, "vs. Controls")) +
        theme(plot.title = element_text(hjust = 0.5, size=22),
                    legend.text = element_text(size = 13),
                    legend.title = element_text(size = 15),
                    axis.title.x = element_text(size = 18), 
                    axis.title.y = element_text(size = 18)) +
        xlab("Log2 Fold Change") +
        ylab("-Log10 p-value") +
        scale_x_continuous(limits = c(log2(xlim_down), log2(xlim_up)))


    
    # Count the number of Up and Down cells for each compartment
    summary <- allcomp %>%
      group_by(compartment, diffranked) %>%
      summarise(cell_count = n(), .groups = "drop") %>%  # Use .groups = "drop" to ungroup after summarise
      spread(key = diffranked, value = cell_count, fill = 0)
    
    # Total Up and Down cells across the entire dataset
    total_up <- sum(summary$Up, na.rm = TRUE)
    total_down <- sum(summary$Down, na.rm = TRUE)
    
    # Proportion of each compartment in the Up and Down cases
    summary <- summary %>%
      group_by(compartment) %>%
      mutate(proportion_up = ifelse(!is.na(Up), Up / total_up, 0),
             proportion_down = ifelse(!is.na(Down), Down / total_down, 0)) %>%
      ungroup()
    
    # Format the output into a nice dataframe
    formatted_summary <- summary %>%
      select(compartment, Up, Down, proportion_up, proportion_down) %>%
      arrange(compartment)
    
    # Add a summary row that sums up the values
    summary_row <- formatted_summary %>%
      summarise(compartment = "Total",
                Up = sum(Up, na.rm = TRUE),
                Down = sum(Down, na.rm = TRUE),
                proportion_up = sum(proportion_up, na.rm = TRUE),
                proportion_down = sum(proportion_down, na.rm = TRUE), .groups = "drop")  # Ungroup the summary row
    
    # Combine the formatted summary with the summary row
    formatted_summary <- bind_rows(formatted_summary, summary_row)
    
    ggsave(sprintf("/plots/volcano_plots/volcano_%s_rna_ct_tissue.pdf", gsub(" ", "_", cancer_type)), plot = plot, width = 11, height = 8)
    return(list(allcomp = allcomp, plot = plot, summary = formatted_summary))
}

result <- plotDifferentialCellTypes(crc_test, "Colorectal Cancer", outdir)

# Set the figure size for R plots in Jupyter
options(repr.plot.width = 10, repr.plot.height = 8)
result[[2]]

volcano_plot_crc <- result[[2]]
volcano_plot_crc

# ggsave("/plots/volcano_crc_rna_ct_tissue.pdf", plot = volcano_plot_crc, width = 11, height = 8)

# Count the number of Up and Down cells for each compartment
summary <- result[[1]] %>%
  group_by(compartment, diffranked) %>%
  summarise(cell_count = n(), .groups = "drop") %>%  # Use .groups = "drop" to ungroup after summarise
  spread(key = diffranked, value = cell_count, fill = 0)

# Total Up and Down cells across the entire dataset
total_up <- sum(summary$Up, na.rm = TRUE)
total_down <- sum(summary$Down, na.rm = TRUE)

# Proportion of each compartment in the Up and Down cases
summary <- summary %>%
  group_by(compartment) %>%
  mutate(proportion_up = ifelse(!is.na(Up), Up / total_up, 0),
         proportion_down = ifelse(!is.na(Down), Down / total_down, 0)) %>%
  ungroup()

# Format the output into a nice dataframe
formatted_summary <- summary %>%
  select(compartment, Up, Down, proportion_up, proportion_down) %>%
  arrange(compartment)

# Add a summary row that sums up the values
summary_row <- formatted_summary %>%
  summarise(compartment = "Total",
            Up = sum(Up, na.rm = TRUE),
            Down = sum(Down, na.rm = TRUE),
            proportion_up = sum(proportion_up, na.rm = TRUE),
            proportion_down = sum(proportion_down, na.rm = TRUE), .groups = "drop")  # Ungroup the summary row

# Combine the formatted summary with the summary row
formatted_summary <- bind_rows(formatted_summary, summary_row)

cancer_types <- list(unique(data$status))[[1]]
cancer_types <- cancer_types[!grepl("Healthy", cancer_types)]

## all cancer types
cancer_types <- list(unique(data$status))[[1]]
cancer_types <- cancer_types[!grepl("Healthy", cancer_types)]

summary_list <- list()

for (cancer_type in cancer_types){
    
    gender_cancer <- data %>% filter(status %in% c(cancer_type)) %>% distinct(gender) %>% pull(gender)
    sub <- data %>% filter(status %in% c(cancer_type, "Healthy") & gender %in% gender_cancer)
    print(cancer_type)
    # print(gender_cancer)
    #print(length(unique(sub$sample)))
    
    table <- plotDifferentialCellTypes(sub, cancer_type, outdir)[1]
    plot <- plotDifferentialCellTypes(sub, cancer_type, outdir)[2]
    summary <- plotDifferentialCellTypes(sub, cancer_type, outdir)[3]
    summary_list[[cancer_type]] <- summary
    # print(summary)
    print(plot)
}

library(viridis)
custom_colors <- viridis(100, option = "viridis")[1:100]
library(patchwork)

for (cancer_type in names(summary_list)){
    down <- summary_list[[cancer_type]][[1]] %>% select(compartment, proportion_down)
    up <- summary_list[[cancer_type]][[1]] %>% select(compartment, proportion_up)
    down$CancerType <- cancer_type
    up$CancerType <- cancer_type
}

# Initialize empty lists to store down and up dataframes
down_list <- list()
up_list <- list()

# Loop through each cancer type to extract and store down and up proportions
for (cancer_type in names(summary_list)) {
  # Extract the down and up proportions for the current cancer type
  down <- summary_list[[cancer_type]][[1]] %>% select(compartment, proportion_down) %>%
    filter(compartment != "Total")  # Remove the 'Total' compartment rows
  
  up <- summary_list[[cancer_type]][[1]] %>% select(compartment, proportion_up) %>%
    filter(compartment != "Total")  # Remove the 'Total' compartment rows
  
  # Add the cancer type as a column
  down$CancerType <- cancer_type
  up$CancerType <- cancer_type
  
  # Append to the list
  down_list[[cancer_type]] <- down
  up_list[[cancer_type]] <- up
}

# Combine all the down and up dataframes into one dataframe each
down_df <- bind_rows(down_list)
up_df <- bind_rows(up_list)

# Initialize empty lists to store down and up dataframes
down_list <- list()
up_list <- list()

# Loop through each cancer type to extract and store down and up proportions
for (cancer_type in names(summary_list)) {
  # Extract the down and up proportions for the current cancer type
  down <- summary_list[[cancer_type]][[1]] %>% select(compartment, proportion_down) %>%
    filter(compartment != "Total")  # Remove the 'Total' compartment rows
  
  up <- summary_list[[cancer_type]][[1]] %>% select(compartment, proportion_up) %>%
    filter(compartment != "Total")  # Remove the 'Total' compartment rows
  
  # Add the cancer type as a column
  down$CancerType <- cancer_type
  up$CancerType <- cancer_type
  
  # Append to the list
  down_list[[cancer_type]] <- down
  up_list[[cancer_type]] <- up
}

# Combine all the down and up dataframes into one dataframe each
down_df <- bind_rows(down_list)
up_df <- bind_rows(up_list)


# Set factor levels for compartments to ensure the correct order
compartment_order <- c("stromal", "endothelial", "epithelial", "immune")

up_down <- list(up_df, down_df)
props <- c("proportion_up", "proportion_down")
titles <- c("Upregulated Cells", "Downregulated Cells")
plots <- list()

compartment_colors <- c(
  "endothelial" = "#FFA500",  # Orange
  "epithelial"  = "#008000",  # Green
  "immune"      = "#0000FF",  # Blue
  "stromal"     = "#999999"   # Grey
)

# Set compartment factor levels for both datasets
up_df$compartment <- factor(up_df$compartment, levels = compartment_order)
down_df$compartment <- factor(down_df$compartment, levels = compartment_order)

for (i in 1:2){
    df <- up_down[[i]]
    prop <- props[i]
    # Set y label conditionally
    y_label <- ifelse(i == 1, "Relative Contribution", "")

    p <- ggplot(df, aes(x = CancerType, y = .data[[prop]], fill = compartment)) +
        geom_bar(stat = "identity", position = "fill") +  
        theme_minimal() +
        scale_fill_manual(values = compartment_colors) +
        labs(
            title = titles[i],
            x = "Cancer Type",
            y = y_label,
            fill = "Compartment"
        ) +
        theme(
            plot.title = element_text(size = 20, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
            axis.text.y = element_text(size = 14),
            axis.title.y = element_text(size = 16),
            axis.title.x = element_blank(),
            legend.position = ifelse(i == 1, "none", "right"),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 14)
        )
        
    plots[[i]] <- p
}

# Combine and display both plots in one figure
combined_plot_compartments <- plots[[1]] + plots[[2]] + plot_layout(ncol = 2)
print(combined_plot_compartments)


# Save the combined plot as a PDF
# ggsave("/plots/volcano_summary_compartments.pdf", combined_plot_compartments, width = 14, height = 6, device = "pdf")

sex_df <- data %>% filter(status == "Healthy") %>% filter(grepl("Uterus", cell_type_tissue))

plotDifferentialCellTypes_sex <- function(data) {
    
    allcomp <- c()
    cell <- unique(data$cell_type_tissue)
    
    #run two-sided wilcox test and compute mean foldchange for each cell type
    for(i in 1:length(cell)) {
    df <- subset(data, cell_type_tissue %in% cell[i])
    pheno.counts <- as.data.table(table(df$gender))
    labels <- paste0(pheno.counts[1,1],".n.", pheno.counts[1,2], ".vrs.", pheno.counts[2,1], ".n.", pheno.counts[2,2])
    df.wilcox <- wilcox_test(df, rank ~ gender, ref.group = pheno.counts[2,1]$V1)
    df.mean <- df %>% group_by(gender) %>% summarize(Mean = mean(rank, na.rm=TRUE)) %>% as.data.frame()
    df.mean.FC <- df.mean$Mean[2]/df.mean$Mean[1] %>% as.data.frame()
    colnames(df.mean.FC) <- c("foldchange")
    df.mean.FC$cell_type_tissue <- cell[i]
    df.all <- cbind(df.wilcox, df.mean.FC)
    allcomp <- rbind(allcomp, df.all)
    }
    
    #multiple testing correction
    allcomp$p.adj <- p.adjust(allcomp$p, method="fdr")

    #order by p.adj
    allcomp <- allcomp[order(allcomp[,10]),]
    
    #create color for volcano plot
    allcomp$diffranked <- "Not Sig"
    allcomp$diffranked[allcomp$foldchange > 1 & allcomp$p.adj < 0.05] <- "Up"
    allcomp$diffranked[allcomp$foldchange < 1 & allcomp$p.adj < 0.05] <- "Down"

    allcomp$delabel <- NA
    allcomp[allcomp$diffranked == "Down" | allcomp$diffranked == "Up",]$delabel <- allcomp[allcomp$diffranked == "Down" | allcomp$diffranked == "Up",]$cell_type

    # add compartment and cell_type_tissue
    allcomp$cell_type_tissue <- str_trim(allcomp$cell_type_tissue)
    allcomp <- allcomp %>% left_join(compartment_map, by = c("cell_type_tissue" = "cell_type_tissue"))
    allcomp <- allcomp %>%
    # Create a new column 'top_label' based on the conditions
    mutate(
    top_label = case_when(
      diffranked == "Down" & foldchange %in% head(sort(foldchange[diffranked == "Down"], decreasing = FALSE), 10) ~ cell_type_tissue,  # Top 10 lowest foldchange for "Down"
      diffranked == "Up" & foldchange %in% head(sort(foldchange[diffranked == "Up"], decreasing = TRUE), 10) ~ cell_type_tissue,    # Top 10 highest foldchange for "Up"
      TRUE ~ NA_character_  # Set NA for all other rows
    )
    )
    return(list(allcomp = allcomp))
}

result <- plotDifferentialCellTypes_sex(sex_df)
sex_healthy <- result[[1]]

plots_sex = list()
for (sex_tissue in c("Prostate", "Uterus")){
    sex_df <- data %>% filter(status == "Healthy") %>% filter(grepl(sex_tissue, cell_type_tissue))
    result <- plotDifferentialCellTypes_sex(sex_df)
    sex_healthy <- result[[1]]

    # Add logFoldChange
    sex_healthy_plot <- sex_healthy %>%
      mutate(logFoldChange = log2(foldchange))
    sex_healthy_plot$cell_type_tissue <- str_remove(sex_healthy_plot$cell_type_tissue, paste0("_", sex_tissue))
    
    custom_colors <- c("Not Sig" = "grey80", "Down" = "#541352FF", "Up" = "#10a53dFF")
    
    
    # Plot with custom gradient
    plot <- ggplot(sex_healthy_plot, aes(x = logFoldChange, y = reorder(cell_type_tissue, logFoldChange), fill = diffranked)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
        theme_minimal() +
        scale_fill_manual(values = custom_colors) +  # Use custom colors
        labs(
            title = paste0(sex_tissue, ": Female vs. Male"),
            x = "Log2 Fold Change",
            y = "", # Cell Type
            fill = "diffranked"
        ) +
        scale_x_continuous(limits = c(-0.15, 0.175)) + 
        theme(
            axis.text.x = element_text(size = 18),  # Increase x-axis labels
            axis.text.y = element_text(size = 18),  # Increase y-axis labels
            axis.title = element_text(size = 18),  # Increase axis titles
            plot.title = element_text(size = 25, hjust = 0.5),  # Increase title size
            legend.text = element_text(size = 16),  # Increase legend text
            legend.title = element_text(size = 18)  # Increase legend title
        )
    print(plot)
    plots_sex[[sex_tissue]] <- plot
}


combined_plot <- plots_sex[["Prostate"]] + plots_sex[["Uterus"]] + plot_layout(nrow = 2, ncol = 1, heights = c(1, 1))

# Display the combined plot
print(combined_plot)


# Save the combined plot as a PDF
# ggsave("/plots/sex_comparison_plot.pdf", combined_plot, width = 13, height = 13, dpi = 300)
