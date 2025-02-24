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
