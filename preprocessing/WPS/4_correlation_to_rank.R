# load libraries
library(Matrix)
library(Seurat)
library(SeuratObject)
library(gplots)

# Define samples for comparison
# samples <- c("EE87891","EE87892","EE87893")
folders = "/mnt/DATA3/timo/main/body"
all_folders <- list.dirs(folders, full.names = TRUE, recursive = FALSE)
samples_to_exclude <- c(
 "EE87920", "EE87921",
"EE87947", "EE87959", "EE87963", "EE87998", 
"EE88009", "EE88015", "EE88020", "EE88025", 
"EE88027", "EE88082", "EE88093", "EE88127", 
"EE88136", "EE88141", "EE88152", "EE88164"
)
samples_raw <- grep("EE", basename(all_folders), value = TRUE)
samples <- setdiff(samples_raw, samples_to_exclude)
print(paste(length(samples), " samples to process."))
# for testing load pre filtered and sampled rds file
# print("Loading .rds test file.")
# tabula_sapiens <- readRDS("/mnt/DATA3/timo/data/sngl_cell_Tabula_Sapiens/seurat_subset_10k_test.rds")
print("Loading Tabula Sapiens.")
tabula_sapiens <- readRDS("/mnt/DATA3/timo/data/sngl_cell_Tabula_Sapiens/Tabula_Sapiens_filtered.rds")
print("File loaded.")

##################################
print("Average Expression per Cell Type Tissue")
avg_expr <- AverageExpression(object = tabula_sapiens,
                       assays = "RNA",
                       group.by = "cell_type_tissue",
                       return.seurat = FALSE)
avg_expr_df <- as.data.frame(avg_expr)
print("Only keep Genes with >= 3 non-zero values.")
# Count non-zero entries for each row directly on the sparse matrix
non_zero_counts <- rowSums(avg_expr_df != 0)
# Create a logical vector to keep rows with at least 3 or more non-zero values
threshold <- 3
keep_rows <- non_zero_counts[non_zero_counts > threshold]
rows_names_keep <- names(keep_rows)
# subset the initial dgCMatrix using the row names (ENSG00010123) which contain at least 3 non-zero values
avg_expr_df_filtered <- avg_expr_df[rownames(avg_expr_df) %in% rows_names_keep, ]

##################################
#fftColumns <- 29:52 # 160-222
print("Set frequencies of interest.")
fftColumns <- 29:52 # 160-222
selFreq <- c("193","196","199")

##################################
# Generating Rank comparison plots
# Replace BH01 with sample you are using as reference in the rank comparison
# reference sample is the fft WPS averaged over all 262 healthy samples
reference_sample <- "reference_244"
print("Loading reference WPS data.")
fdata <- read.table(sprintf("/mnt/DATA3/timo/main/body/fft_summaries/fft_%s_WPS.tsv.gz", reference_sample),as.is=T,sep="\t",header=T,comment.char="~")
fdata <- cbind(X.Region = rownames(fdata), fdata)
colnames(fdata) <- sub("X","",colnames(fdata))
rownames(fdata) <- fdata[,1]
fdata <- fdata[,c(1,rev(c(2:dim(fdata)[2])))]
# get common rows
common_rownames <- intersect(rownames(fdata), rownames(avg_expr_df_filtered))
# select common rows
print(paste("Number of ENSGs: ", length(common_rownames)))
fdata_subset <- fdata[common_rownames, , drop = FALSE]
avg_expr_df_filtered_subset <- avg_expr_df_filtered[common_rownames, , drop = FALSE]

print("Compute Correlation ranks.")
refCorrelation <- cor(rowMeans(fdata_subset[,selFreq]),avg_expr_df_filtered_subset[,order(names(avg_expr_df_filtered_subset))],use="pairwise.complete.obs")
print(paste("Reference Correlation for ", reference_sample, "calulated."))

#pdf("Ave193-199bp_correlation_rank.pdf",width=10,height=15)
print("Going through samples...")
for (sample in samples){
    print(paste("Ranking ", sample, "to reference sample."))
    fdata <- read.table(sprintf("/mnt/DATA3/timo/main/body/fft_summaries/fft_%s_WPS.tsv.gz",sample),as.is=T,sep="\t",header=T,comment.char="~")
    colnames(fdata) <- sub("X","",colnames(fdata))
    rownames(fdata) <- fdata[,1]
    fdata <- fdata[,c(1,rev(c(2:dim(fdata)[2])))]
    # get common rows
    common_rownames <- intersect(rownames(fdata), rownames(avg_expr_df_filtered))
    # select common rows
    fdata_subset <- fdata[common_rownames, , drop = FALSE]
    avg_expr_df_filtered_subset <- avg_expr_df_filtered[common_rownames, , drop = FALSE]

    res <- cor(rowMeans(fdata_subset[,selFreq]),avg_expr_df_filtered_subset[,order(names(avg_expr_df_filtered_subset))],use="pairwise.complete.obs")
    #res <- data.frame(category=tLabels$Category,description=tLabels$Type,tissue=colnames(res),correlation=as.numeric(res),rankDiff=rank(refCorrelation)-rank(res))
    res <- data.frame(cell_type=colnames(res),correlation_sample=as.numeric(res),correlation_ref=as.numeric(refCorrelation), rank_sample=rank(res), rank_ref=rank(refCorrelation),rankDiff=rank(refCorrelation)-rank(res))
    res <- res[order(res$rankDiff, decreasing = TRUE), ]
    write.csv(res, sprintf("/mnt/DATA3/timo/main/body/outputs/corr_tissues/244/%s_244_correlation.csv", sample), row.names = FALSE)
    }
print("Done.")
