library(data.table)
library(tidyverse)
library(e1071)
library(pROC)

myseed <- 42 ##Set seed
cost <- 1 ##SVM regularization, cost of constraints

# Load data
data <- read.csv("data/ranks/Ranks_FCC_RNA_cell_type_tissue.csv", header = TRUE)

cancer_types = list(unique(data$status))[[1]]
cancer_types = cancer_types[cancer_types != "Healthy"]

# Set seed for reproducibility
set.seed(myseed)

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

# Loop through each cancer type for One-vs-Rest classification
for (cancer_type in cancer_types) {
    # print(paste("Training OVR SVM for:", cancer_type))
    
    # Filter gender for Ovarian Cancer (only use female samples)
    if (cancer_type %in% c("Ovarian Cancer", "Breast Cancer")) {
        sub <- data %>% filter(gender == "F")  # Subset for female samples
    } else {
        sub <- data  # No gender filter for other cancer types
    }

    # Prepare dataset: Remove unnecessary columns (gender, sample, tfx, correlation)
    combined <- sub
    combined$group <- paste0(combined$sample, ":", combined$status)
    
    # Pivot wider to turn cell type-tissue into columns
    combined.wider <- combined %>%
    select(-sample, -tfx, -gender, -correlation) %>%
    pivot_wider(names_from = cell_type_tissue, values_from = rank) %>%
    as.data.frame()
    
    # Set rownames to group
    rownames(combined.wider) <- combined.wider$group
    
    # Convert labels to binary for One-vs-Rest (1 for current cancer type, 2 for all others)
    dtmat_group <- ifelse(combined.wider$status == cancer_type, 1, 2)
    dtmat_group <- as.factor(dtmat_group)

    dtmat <- combined.wider %>% select(-group, -status) %>% data.matrix()

    # Check class distribution in the training set
    # print(paste("Class distribution for", cancer_type, ":", table(dtmat_group)))
    
    # Class weights for handling imbalance
    wts <- 100 / table(dtmat_group)
    n_train <- nrow(dtmat)
    # print(paste("Class weights for", cancer_type, ":", wts))

    n_features <- ncol(dtmat)
    feature_weights <- matrix(0, nrow = n_features, ncol = n_train)
    rownames(feature_weights) <- colnames(dtmat)
    colnames(feature_weights) <- paste0("LeaveOut_", 1:n_train)
    
    # Initialize variables to store results for LOO-CV
    alltruegroup <- c()
    allpredgroup <- c()
    allpredval <- c()

    # Perform Leave-One-Out Cross Validation
    for (i in 1:n_train) {
        # Define the test set (leave one out)
        test_idx <- i
        train_idx <- setdiff(1:nrow(dtmat), test_idx)
        
        traindt <- dtmat[train_idx, ]
        traingroup <- dtmat_group[train_idx]
        testdt <- dtmat[test_idx, , drop = FALSE]
        testgroup <- dtmat_group[test_idx]
        
        # Train SVM with linear kernel on the training set
        svm.model <- svm(
          traindt, traingroup,
          class.weights = wts, kernel = "linear", cost = cost,
          scale = FALSE, decision.values = TRUE
        )
        
        w <- t(svm.model$coefs) %*% svm.model$SV
        feature_weights[, i] <- as.numeric(w)
        
        svm.pred <- predict(svm.model, testdt, decision.values = TRUE)
        
        # Make predictions on the test set
        # numeric_true <- ifelse(dtmat_group[i] == cancer_type, 1, 2)
        alltruegroup <- c(alltruegroup, dtmat_group[i])
        allpredgroup <- c(allpredgroup, attributes(svm.pred)$decision.values)
        # numeric_pred <- ifelse(svm.pred == cancer_type, 1, 2)
        allpredval <- c(allpredval, svm.pred)
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

    alltruegroup <- factor(alltruegroup, levels = c(2, 1), labels = c("Control", "Cancer Type"))
    allpredval <- factor(allpredval, levels = c(2, 1), labels = c("Control", "Cancer Type"))

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
    plot(svmperfm, print.auc = TRUE, main = paste("ROC for", cancer_type))
    text(x = 0.35, y = 0.1, labels = sprintf(
      "Accuracy:      %.2f\nSensitivity:    %.2f\nSpecificity:    %.2f\nPrecision:      %.2f\nF1 Score:      %.2f",
      round(accuracy, 2), round(sensitivity, 2), round(specificity, 2), round(precision, 2), round(f1_score, 2)),
      pos = 4, cex = 1.0, col = "black")
    
}

feature_weights_df <- bind_rows(feature_weights_dict)

# Extract AUC values from the list
auc_values <- map_dbl(roc_list, ~ .x$auc)

# Reorder roc_list by descending AUC
roc_list_sorted <- roc_list[order(auc_values, decreasing = FALSE)]

# pdf("/mnt/DATA3/timo/plots/roc_curves_OvR.pdf", width = 11, height = 11)

colors <- c("#1F77B4", "#E377C2", "#B8860B", "#FF7F0E", "#D62728", "#2CA02C", "#6A0DAD")

# par(mar = c(10, 10, 1, 2))

# Create a blank plot with the first ROC curve
plot(roc_list_sorted[[1]], col=colors[1], main="ROC Curves of OvR SVM", #"ROC Curves for all Cancer Types"
    xlim=c(1,0), ylim=c(0,1), 
    cex.main=3,
    cex.lab=2.25, # 2
    cex.axis=1.5, # 1.75
    font.main = 1)#,
    # col.lab = "white")

# Manually add custom axes with step size 0.5
# axis(1, at=seq(1, 0, by=-.2), cex.axis=2)  # X-axis (Specificity)
# axis(2, at=seq(0, 1, by=.2), cex.axis=2, las = 1)  # Y-axis (Sensitivity)

# # Add custom axis labels
# mtext("Specificity", side=1, line=5, cex=2.5)  
# mtext("Sensitivity", side=2, line=5, cex=2.5)  
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

metrics_table$ref_data_type = "rna_ct_tissue_ovr"
# Save the table to a file
# write.csv(metrics_table, file = "/outputs/eva_metrics/rna_ct_tissue_eva_metrics_table_OvR.csv", row.names = FALSE)
