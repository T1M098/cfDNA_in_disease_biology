library(data.table)
library(tidyverse)
library(e1071)
library(pROC)
args <- commandArgs(TRUE)

myseed <- 42
cost <- 1 

# Load data
refdata_type <- "RNA" # "ATAC" or "RNA"
data <- read.csv(sprintf("/data/ranks/Ranks_FCC_%s_cell_type.csv",  refdata_type), header = TRUE)
head(data)

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

# Initialize an empty data frame to store metrics
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

for (cancer_type in cancer_types) {
  print(cancer_type)
  gender_cancer <- data %>% filter(status == cancer_type) %>% pull(gender) %>% unique()
  sub <- data %>% filter(
    status %in% c(cancer_type, "Healthy") &
    gender %in% gender_cancer)

  combined <- sub
  combined$group <- paste0(combined$sample, ":", combined$status)
  combined.wider <- combined %>% select(-sample, -status, -tfx, -gender, -correlation) %>% pivot_wider(names_from = cell_type, values_from = rank) %>% as.data.frame()
  rownames(combined.wider) <- combined.wider$group
  dtmat <- combined.wider %>% select(-group) %>% data.matrix()

  # Set seed
  set.seed(myseed)

  # Set variables
  orgroup <- row.names(dtmat)
  dtmat_group <- sapply(strsplit(row.names(dtmat), split=':', fixed=TRUE), function(x) (x[2]))
  dtmat_group <- as.factor(dtmat_group)

  # Get weights
  wts <- 100 / table(dtmat_group)
  n_train <- nrow(dtmat)

  alltruegroup <- c()
  allpredgroup <- c()
  allpredval <- c()
  predres <- data.frame(Group = character(), DecisionValue = numeric())

  # Run leave-one-out SVM
  for (x in 1:n_train) {
    traindt <- dtmat[-x,]
    traingroup <- dtmat_group[-x]
    testdt <- t(dtmat[x,])
    testgroup <- dtmat_group[x]
    svm.model <- svm(traindt, traingroup, class.weights = wts, kernel="linear", cost=cost, scale=FALSE, decision.values=TRUE)
    svm.pred <- predict(svm.model, testdt, decision.values=TRUE)
    numeric_true <- ifelse(testgroup == "Healthy", 2, 1)
    alltruegroup <- c(alltruegroup, numeric_true)
    allpredgroup <- c(allpredgroup, attributes(svm.pred)$decision.values)
    numeric_pred <- ifelse(svm.pred == "Healthy", 2, 1)
    allpredval <- c(allpredval, numeric_pred)
  }

  alltruegroup <- factor(alltruegroup, levels = c(2, 1), labels = c("Control", "Cancer"))
  allpredval <- factor(allpredval, levels = c(2, 1), labels = c("Control", "Cancer"))

  # Calculate performance metrics
  svmperfm <- roc(response=alltruegroup, predictor=allpredgroup, ci=TRUE)
  auc_value <- svmperfm$auc

  runacc <- table(alltruegroup, allpredval)
  TP <- runacc[2,2]
  TN <- runacc[1,1]
  FP <- runacc[1,2]
  FN <- runacc[2,1]

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
  plot(svmperfm, print.auc=TRUE, main = paste("ROC ", cancer_type))
  text(x = 0.35, y = 0.1, labels = sprintf(
    "Accuracy:      %.2f\nSensitivity:    %.2f\nSpecificity:    %.2f\nPrecision:      %.2f\nF1 Score:      %.2f",
    round(accuracy, 2), round(sensitivity, 2), round(specificity, 2), round(precision, 2), round(f1_score, 2)),
    pos = 4, cex = 1.0, col = "black")
}

metrics_table$ref_data_type = refdata_type

# Save the table to a file
# write.csv(metrics_table, file = sprintf("/ATAC/outputs/eva_metrics/%s_eva_metrics_table.csv", refdata_type), row.names = FALSE)
