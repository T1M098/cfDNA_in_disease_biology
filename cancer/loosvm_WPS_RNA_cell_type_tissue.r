library(data.table)
library(tidyverse)
library(e1071)
library(pROC)
args <- commandArgs(TRUE)

myseed <- 42
cost <- 1

# Load data
dat <- read.csv("/data/ranks/Ranks_WPS_RNA_cell_type_tissue.tsv", header = TRUE, sep = "\t")

# Remove non cfDNA samples
noncfDNA_samples <- c("EE87920", "EE87921",
    "EE87947", "EE87959", "EE87963", "EE87998", 
    "EE88009", "EE88015", "EE88020", "EE88025", 
    "EE88027", "EE88082", "EE88093", "EE88127", 
    "EE88136", "EE88141", "EE88152", "EE88164"
    )

dat <- na.omit(dat)

dat <- dat %>% filter(!sample %in% noncfDNA_samples)
# Remove the "RNA." prefix from the cell_type column
dat$cell_type <- gsub("^RNA\\.", "", dat$cell_type)
 
meta_dat <- fread("/data/Cristiano_samplemap.tsv", header = TRUE)
gender <- meta_dat[,c("Gender", "sample")]

# group data according to cancer types
colorectal <- dat[dat$status == "Colorectal Cancer",]
breast <- dat[dat$status == "Breast Cancer",]
lung <- dat[dat$status == "Lung Cancer",]
ovarian <- dat[dat$status == "Ovarian Cancer",]
pancreatic <- dat[dat$status == "Pancreatic Cancer",]
duedeonal <- dat[dat$status == "Duodenal Cancer",]
bile <- dat[dat$status == "Bile Duct Cancer",]
gastric <- dat[dat$status == "Gastric cancer",]

# also group match healthy cases i.e. breast cancer, ovarian cancer only with healthy female samples
healthy_females <- dat %>% 
    filter(status == "Healthy") %>%
    inner_join(gender %>% filter(Gender == "F"), by = "sample") %>% select(-Gender)
all_healthy <- dat[dat$status == "Healthy",]

healthy_sample_ids <- unique(all_healthy$sample)
healthy_female_sample_ids <- unique(healthy_females$sample)

set.seed(42)
# add healthy samples
sample_ids_healthy <- sample(healthy_sample_ids, 20)
sample_ids_healthy_females <- sample(healthy_female_sample_ids, 20)

colorectal <- rbind(colorectal, all_healthy)# %>% filter(sample %in% sample_ids_healthy))
lung <- rbind(lung, all_healthy)
pancreatic <- rbind(pancreatic, all_healthy)
duedeonal <- rbind(duedeonal, all_healthy)
bile <- rbind(bile, all_healthy)
gastric <- rbind(gastric, all_healthy)
# female
ovarian <- rbind(ovarian, healthy_females)
breast <- rbind(breast, healthy_females)

dic_ct <- list(
  'colorectal' = colorectal,
  'breast' = breast,
  'lung' = lung,
  'ovarian' = ovarian,
  'pancreatic' = pancreatic,
  'bile' = bile,
  'gastric' = gastric
)

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

# Predictions for all cancer types:
for (cancer_type in names(dic_ct)){
    print(cancer_type)

    # testing with normalized correlations values instead of ranks as features
    combined <- dic_ct[[cancer_type]]
    combined$group <- paste0(combined$sample, ":", combined$status)
    combined.wider <- combined %>% select(-sample, -status) %>% pivot_wider(names_from = cell_type, values_from = rank_sample) %>% as.data.frame()
    rownames(combined.wider) <- combined.wider$group
    dtmat <- combined.wider %>% select(-group) %>% data.matrix()

    #set seed
    set.seed(myseed)

    #set variables
    orgroup <- row.names(dtmat)
    dtmat_group <- sapply(strsplit(row.names(dtmat), split=':', fixed=TRUE), function(x) (x[2]))
    dtmat_group <- as.factor(dtmat_group)

    #get weights
    wts <- 100 / table(dtmat_group)
    n_train <- nrow(dtmat)

    alltruegroup <- c()
    allpredgroup <- c()
    allpredval <- c()
    predres <- c()

    #run leave-one-out svm
    for (x in 1:n_train) {
      traindt <- dtmat[-x,]
      traingroup <- dtmat_group[-x]
      testdt <- t(dtmat[x,])
      #train.pca <- prcomp(traindt)
      #print(dim(train.pca$x))
      #project.test <- predict(train.pca, testdt)
      testgroup <- dtmat_group[x]
      svm.model <- svm(traindt, traingroup, class.weights = wts, kernel="linear", cost=cost, scale=FALSE, decision.values=TRUE)
      svm.pred <- predict(svm.model,testdt,decision.values=TRUE)
      numeric_true <- ifelse(testgroup == "Healthy", 2, 1)
      alltruegroup <- c(alltruegroup,numeric_true)
      allpredgroup <- c(allpredgroup,attributes(svm.pred)$decision.values)
      numeric_pred <- ifelse(svm.pred == "Healthy", 2, 1)
      allpredval <- c(allpredval,numeric_pred)
      predres <- rbind(predres, data.frame(orgroup[x],attributes(svm.pred)$decision.values))
      #print(data.frame(orgroup[x],svm.pred))
    }
            
    # Assigne control and case levels to 1 and 2
    alltruegroup <- factor(alltruegroup, levels = c(2, 1), labels = c("Control", "Cancer"))
    allpredval <- factor(allpredval, levels = c(2, 1), labels = c("Control", "Cancer"))
    #export decision values and print performance metrics
    #write.table(predres, paste0(output_dir, "/svm_loo_decision_values.txt", sep="\t", row.names=F, col.names=T, quote=F)
    svmperfm <- roc(response=alltruegroup,predictor=allpredgroup,ci=TRUE)
    auc_value <- svmperfm$auc
    print(svmperfm$auc)
    print(svmperfm$ci)
    runacc <- table(alltruegroup,allpredval)
    print(runacc)

    TP <- runacc[2,2]
    TN <- runacc[1,1]
    FP <- runacc[1,2]
    FN <-runacc[2,1]
                          
    # Calculate Accuracy
    accuracy <- (TP + TN) / (TP + TN + FP + FN)

    # Calculate Sensitivity (Recall)
    sensitivity <- TP / (TP + FN)

    # Calculate Specificity
    specificity <- TN / (TN + FP)

    # Calculate Precision
    precision <- TP / (TP + FP)

    # Calculate F1 Score
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

    #plot ROC curve
    plot(svmperfm,print.auc=TRUE, main = paste("ROC ", cancer_type))
    # Add performance metrics to the plot
    text(x = 0.35, y = 0.1, labels = sprintf(
        "Accuracy:      %.2f\nSensitivity:    %.2f\nSpecificity:    %.2f\nPrecision:      %.2f\nF1 Score:      %.2f",
        round(accuracy,2), round(sensitivity,2), round(specificity, 2), round(precision, 2), round(f1_score, 2)),
        pos = 4, cex = 1.0, col = "black")
    #dev.off()
}

metrics_table$ref_data_type = "WPS RNA Cell Type Tissue"
metrics_table <- metrics_table %>%
  mutate(CancerType = case_when(
    CancerType == "colorectal" ~ "Colorectal Cancer",
    CancerType == "breast" ~ "Breast Cancer",
    CancerType == "lung" ~ "Lung Cancer",
    CancerType == "ovarian" ~ "Ovarian Cancer",
    CancerType == "pancreatic" ~ "Pancreatic Cancer",
    CancerType == "bile" ~ "Bile Duct Cancer",
    CancerType == "gastric" ~ "Gastric cancer",
    TRUE ~ CancerType  # Keep other values unchanged
  ))

# Save the table to a file
# write.csv(metrics_table, file = "/outputs/eva_metrics/WPS_RNA_cell_type_tissue_eva_metrics_table.csv", row.names = FALSE)
