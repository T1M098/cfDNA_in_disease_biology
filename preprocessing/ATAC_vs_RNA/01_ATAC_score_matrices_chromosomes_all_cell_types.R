library(dplyr)
library(data.table)
library(readr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
chromosome_input <- args[1]

print(sprintf("Compute ATAC score matrix for chromosome %s.",chromosome_input))


# Load gene mapping
print("Loading Gene Map.")
gene_mapping <- read.table("/mnt/DATA3/timo/data/annotation_files/final_annotation_file.tsv",
                           sep = "\t", header = FALSE, 
                           col.names = c("GeneID", "chromosome", "start", "end", "strand"))


# Filter out unwanted chromosomes
gene_mapping_cleaned <- gene_mapping %>%
  filter(!chromosome %in% c("M", "X", "Y"))

# Drop duplicate GeneIDs, keeping the first occurrence
gene_mapping_cleaned <- gene_mapping_cleaned %>%
  distinct(GeneID, .keep_all = TRUE)

# Convert 'start' and 'end' columns to numeric
gene_mapping_cleaned <- gene_mapping_cleaned %>%
  mutate(across(c(start, end), as.numeric))

# Filter annotation file, chromosome == chromosome_input
gene_mapping_cleaned <- gene_mapping_cleaned %>%
  filter(chromosome == as.character(chromosome_input))

print("Loading ATAC Data.")
atac <- fread("/mnt/DATA1/Zhang/yv4fzv6cnm-4/Figure_2/2E_Accessibility_score.tsv")

df_atac <- as.data.frame(atac)
    
atac_subset <- df_atac

# Check the dimensions
dim(atac_subset)

chrm_string <- paste0("chr", chromosome_input, ":")

# Filter rows containing desired chromosome
atac_subset <- atac_subset %>%
  filter(grepl(chrm_string, V1))

atac_subset <- atac_subset %>%
    separate(V1, into = c("chromosome", "start_end"), sep = ":") %>%
    separate(start_end, into = c("start", "end"), sep = "-")

atac_subset$start <- as.numeric(atac_subset$start)
    
atac_subset$end <- as.numeric(atac_subset$end)
    atac_subset$chromosome <- gsub("^chr", "", atac_subset$chromosome)


print("Dimension:")
print(dim(atac_subset))
      
atac_subset <- atac_subset %>%
  mutate(GeneID = NA) %>%
  relocate(GeneID, .before = chromosome)
    
atac_subset <- atac_subset %>%
  mutate(start_mapping = NA,
         end_mapping = NA,
         overlap = NA)
    
atac_subset <- atac_subset %>%
    relocate(start_mapping, .after = end) %>%
    relocate(end_mapping, .after = start_mapping) %>%
    relocate(overlap, .after = end_mapping) 

bin_df <- data.frame(matrix(ncol = length(colnames(atac_subset)), nrow = 0))
    
colnames(bin_df) <- colnames(atac_subset)

bin_size <- 400

print("Map ATAC regions to Genes.")
for (i in 1:nrow(atac_subset)){
    start <- atac_subset$start[i]
    end <- atac_subset$end[i]
    overlap_gene <- gene_mapping_cleaned$GeneID[
        (start >= gene_mapping_cleaned$start & end <= gene_mapping_cleaned$end) |  # Fully contained
        (start <= gene_mapping_cleaned$end & end >= gene_mapping_cleaned$start)]    # Partial overlap (either start or end overlaps)
    
    for (gene in overlap_gene){
        atac_subset$GeneID[i] <- gene
        atac_subset$start_mapping[i] <- gene_mapping_cleaned$start[gene_mapping_cleaned$GeneID == gene]
        atac_subset$end_mapping[i] <- gene_mapping_cleaned$end[gene_mapping_cleaned$GeneID == gene]
        # Length of overlap
        overlap_start <- pmax(start, atac_subset$start_mapping[i])  
        overlap_end <- pmin(end, atac_subset$end_mapping[i])
        # Proportion of overlap
        atac_subset$overlap[i] <-  (overlap_end - overlap_start)/bin_size
        bin_df <- rbind(bin_df, atac_subset[i,])
    }

}

# Multiply all accessibility score with the proportion of the overlap (usually it's 1)
bin_df <- bin_df %>%
  mutate(across((which(names(bin_df) == "overlap") + 1):ncol(bin_df), ~ . * overlap))

# Summing columns from the 7th column onward, grouped by GeneID
print("Summarize genewise.")
df_summarized <- bin_df %>%
  group_by(GeneID) %>%
  summarise(across((which(names(bin_df) == "overlap")):(ncol(bin_df) - 1), 
                   \(x) sum(x, na.rm = TRUE)))
final_df <- df_summarized %>% arrange(GeneID)

print("Save outputfile.")
write.table(final_df, sprintf("/mnt/DATA3/timo/ATAC/atac_score_matrices_all_cell_types/atac_scores_matrix_chr_%s.tsv", chromosome_input), sep = "\t", row.names = FALSE, quote = FALSE)


