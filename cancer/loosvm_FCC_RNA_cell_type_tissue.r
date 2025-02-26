library(data.table)
library(tidyverse)
library(e1071)
library(pROC)
args <- commandArgs(TRUE)

myseed <- 42
cost <- 1

# Load data
data <- read.csv("/data/ranks/Ranks_FCC_RNA_cell_type_tissue.csv", header = TRUE)

cancer_types = list(unique(data$status))[[1]]
cancer_types = cancer_types[cancer_types != "Healthy"]

cancer_type = cancer_types[2]
gender_cancer = data %>% filter(status == cancer_type) %>% pull(gender) %>% unique()
sub <- data %>% filter(
    status %in% c(cancer_type, "Healthy") &
    gender %in% gender_cancer)


gender_cancer = data %>% filter(status == cancer_type) %>% pull(gender) %>% unique()
sub <- data %>% filter(
status %in% c(cancer_type, "Healthy") &
gender %in% gender_cancer)

# Initialize an empty metrics table
metrics_table <- data.frame(
  CancerType = character(),
  AUC = numeric(),
  Accuracy = numeric(),
  Sensitivity = numeric(),
  Specificity = numeric(),
  Precision = numeric(),
  F1_Score = numeric(),
  stringsAsFactors = FALSE
)

# Create an empty list to store the ROC objects for each cancer type
roc_list <- list()

# Initialize an empty dictionary to store feature weights
feature_weights_dict <- list()

for (cancer_type in cancer_types) {
    print(cancer_type)
    
    gender_cancer <- data %>% filter(status == cancer_type) %>% pull(gender) %>% unique()
    sub <- data %>% filter(
        status %in% c(cancer_type, "Healthy") &
        gender %in% gender_cancer
        )

    combined <- sub
    combined$group <- paste0(combined$sample, ":", combined$status)
    combined.wider <- combined %>%
        select(-sample, -status, -tfx, -gender, -correlation) %>%
        pivot_wider(names_from = cell_type_tissue, values_from = rank, values_fill = list(rank = 0)) %>%
        as.data.frame()
    
    rownames(combined.wider) <- combined.wider$group
    dtmat <- combined.wider %>%
        select(-group) %>%
        data.matrix()

    set.seed(myseed)
    
    orgroup <- row.names(dtmat)
    dtmat_group <- sapply(strsplit(row.names(dtmat), split = ":", fixed = TRUE), function(x) (x[2]))
    dtmat_group <- as.factor(dtmat_group)
    
    wts <- 100 / table(dtmat_group)
    n_train <- nrow(dtmat)
    
    n_features <- ncol(dtmat)
    feature_weights <- matrix(0, nrow = n_features, ncol = n_train)
    rownames(feature_weights) <- colnames(dtmat)
    colnames(feature_weights) <- paste0("LeaveOut_", 1:n_train)
    
    alltruegroup <- c()
    allpredgroup <- c()
    allpredval <- c()
  
  for (x in 1:n_train) {
    traindt <- dtmat[-x, ]
    traingroup <- dtmat_group[-x]
    testdt <- t(as.matrix(dtmat[x, ]))
    
    if (length(unique(traingroup)) < 2) {
      warning("Not enough groups for classification. Skipping this iteration.")
      next
    }
    
    tryCatch({
      svm.model <- svm(
        traindt, traingroup, 
        class.weights = wts, kernel = "linear", cost = cost, 
        scale = FALSE, decision.values = TRUE
      )
      
      w <- t(svm.model$coefs) %*% svm.model$SV
      feature_weights[, x] <- as.numeric(w)
      
      svm.pred <- predict(svm.model, testdt, decision.values = TRUE)
      
      numeric_true <- ifelse(dtmat_group[x] == "Healthy", 2, 1)
      alltruegroup <- c(alltruegroup, numeric_true)
      allpredgroup <- c(allpredgroup, attributes(svm.pred)$decision.values)
      numeric_pred <- ifelse(svm.pred == "Healthy", 2, 1)
      allpredval <- c(allpredval, numeric_pred)
      
    }, error = function(e) {
      warning("SVM failed for iteration ", x, ": ", e$message)
    })
  }
  
    avg_weights <- rowMeans(feature_weights, na.rm = TRUE)
        important_features <- data.frame(
        Feature = colnames(dtmat),
        AvgWeight = avg_weights
        ) %>% arrange(desc(AvgWeight))
  
    # Add the cancer type to the dataframe
    important_features$CancerType <- cancer_type
    
    # Save to the dictionary
    feature_weights_dict[[cancer_type]] <- important_features
    
    alltruegroup <- factor(alltruegroup, levels = c(2, 1), labels = c("Control", "Cancer"))
    allpredval <- factor(allpredval, levels = c(2, 1), labels = c("Control", "Cancer"))
    
    # Calculate performance metrics
    svmperfm <- roc(response = alltruegroup, predictor = allpredgroup, ci = TRUE)
    auc_value <- svmperfm$auc
    
    runacc <- table(alltruegroup, allpredval)
    print(runacc)
    TP <- runacc[2, 2]
    TN <- runacc[1, 1]
    FP <- runacc[1, 2]
    FN <- runacc[2, 1]
  
    accuracy <- (TP + TN) / (TP + TN + FP + FN)
    sensitivity <- TP / (TP + FN)
    specificity <- TN / (TN + FP)
    precision <- TP / (TP + FP)
    f1_score <- 2 * (precision * sensitivity) / (precision + sensitivity)
    
    print(paste0("Accuracy: ", accuracy))
    
    # Store metrics in the table
    metrics_table <- rbind(metrics_table, data.frame(
        CancerType = cancer_type,
        AUC = auc_value,
        Accuracy = accuracy,
        Sensitivity = sensitivity,
        Specificity = specificity,
        Precision = precision,
        F1_Score = f1_score
    ))

    # Add the ROC curve to the list
    roc_list[[cancer_type]] <- svmperfm
  
    # Plot ROC curve
    plot(svmperfm, print.auc = TRUE, main = paste("ROC ", cancer_type))
        text(x = 0.35, y = 0.1, labels = sprintf(
        "Accuracy:      %.2f\nSensitivity:    %.2f\nSpecificity:    %.2f\nPrecision:      %.2f\nF1 Score:      %.2f",
        round(accuracy, 2), round(sensitivity, 2), round(specificity, 2), round(precision, 2), round(f1_score, 2)),
        pos = 4, cex = 1.0, col = "black")
}


metrics_table$ref_data_type = "rna_ct_tissue"

# Save the table to a file
# write.csv(metrics_table, file = "/outputs/eva_metrics/FCC_RNA_cell_type_tissue_eva_metrics_table.csv", row.names = FALSE)

# Extract AUC values from the list
auc_values <- map_dbl(roc_list, ~ .x$auc)

# Reorder roc_list by descending AUC
roc_list_sorted <- roc_list[order(auc_values, decreasing = FALSE)]

# pdf("/mnt/DATA3/timo/plots/roc_curves_all_cancer_types.pdf", width = 11, height = 11)

colors <- c("#E377C2", "#B8860B", "#2CA02C", "#1F77B4", "#D62728", "#FF7F0E", "#6A0DAD")

# par(mar = c(10, 10, 1, 2))

# Create a blank plot with the first ROC curve
plot(roc_list_sorted[[1]], col=colors[1], main="ROC Curves for all Cancer Types", #"ROC Curves for all Cancer Types"
    xlim=c(1,0), ylim=c(0,1), 
    cex.main=3,
    cex.lab=2.25, # 2
    cex.axis=1.5, # 1.75
    font.main = 1)#,
    # col.lab = "white")

# Loop through all the cancer types and add ROC curves with AUC labels
for (i in 1:length(roc_list_sorted)) {
  # Add ROC curve to the plot
  plot(roc_list_sorted[[i]], col=colors[i], add=TRUE, lwd = 5)
  
  # Calculate the AUC value and its confidence interval
  auc_value <- round(roc_list_sorted[[i]]$auc, 3)
  auc_ci_lower <- round(roc_list_sorted[[i]]$ci[1], 3)
  auc_ci_upper <- round(roc_list_sorted[[i]]$ci[3], 3)
  
  # Place the AUC and CI values on the right-hand bottom
  text(x = 0.8, y = 0.05 + (i - 1) * 0.05, 
       labels = paste0(names(roc_list_sorted)[i], 
                      " AUC: ", auc_value, 
                      " (", auc_ci_lower, " - ", auc_ci_upper, ")"), 
       col=colors[i], cex=2, pos=4)
}

# Optional: Add an axis and gridlines for better readability
abline(h=0, v=0, col="gray", lty=2)

# dev.off()
