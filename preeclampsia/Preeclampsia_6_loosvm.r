library(data.table)
library(tidyverse)
library(e1071)
library(pROC)
args <- commandArgs(TRUE)

myseed <- 42
cost <- 1

# Load data
data <- read.csv("/data/ranks/Ranks_FCC_RNA_cell_type_tissue_preeclampsia.csv", header = TRUE)
data <- na.omit(data)

combined <- data

combined$group <- paste0(combined$sample, ":", combined$status)
combined.wider <- combined %>% select(-sample, -status, -correlation, -ALT) %>% pivot_wider(names_from = cell_type_tissue, values_from = rank) %>% as.data.frame()
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

feature_weights_dict <- list()
n_features <- ncol(dtmat)
feature_weights <- matrix(0, nrow = n_features, ncol = n_train)
rownames(feature_weights) <- colnames(dtmat)
colnames(feature_weights) <- paste0("LeaveOut_", 1:n_train)

#n_train <- 50                   
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
    
    w <- t(svm.model$coefs) %*% svm.model$SV
    feature_weights[, x] <- as.numeric(w)
    
    svm.pred <- predict(svm.model,testdt,decision.values=TRUE)
    numeric_true <- ifelse(testgroup == "Healthy", 2, 1)
    alltruegroup <- c(alltruegroup,numeric_true)
    allpredgroup <- c(allpredgroup,attributes(svm.pred)$decision.values)
    numeric_pred <- ifelse(svm.pred == "Healthy", 2, 1)
    allpredval <- c(allpredval,numeric_pred)
    predres <- rbind(predres, data.frame(orgroup[x],attributes(svm.pred)$decision.values))
    #print(data.frame(orgroup[x],svm.pred))
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

# Assigne control and case levels to 1 and 2
alltruegroup <- factor(alltruegroup, levels = c(2, 1), labels = c("Control", "Preeclampsia"))
allpredval <- factor(allpredval, levels = c(2, 1), labels = c("Control", "Preeclampsia"))


#export decision values and print performance metrics
# Convert true groups to factors with specified levels
#alltruegroup <- factor(alltruegroup, levels = c("Cancer", "Control"))  # Make sure the correct labels are used
#write.table(predres, paste0(output_dir, "/svm_loo_decision_values.txt", sep="\t", row.names=F, col.names=T, quote=F)
svmperfm <- roc(response=alltruegroup,predictor=allpredgroup,ci=TRUE)
print(svmperfm$auc)
print(svmperfm$ci)
runacc <- table(alltruegroup,allpredval)
print(runacc)
acclist <- (runacc[1,1]+runacc[2,2])/sum(runacc)
print(paste0("Accuracy: ", acclist))

#plot ROC curve
#pdf(paste0(output_dir, "/svm_loo_roc.pdf"))
plot(svmperfm,print.auc=TRUE)
#dev.off()

feature_weights_df <- bind_rows(feature_weights_dict)

# Save as CSV
# write.csv(feature_weights_df, "/data/feature_weights/FCC_RNA_cell_type_tissue_loosvm_feature_weights_preeclampsia.csv", row.names = FALSE)

preeclampsia <- data
dic_ct <- list(
  'preeclampsia' = data
)

# Predictions for all cancer types:
for (cancer_type in names(dic_ct)){
    print(cancer_type)

    # testing with normalized correlations values instead of ranks as features
    combined <- dic_ct[[cancer_type]]
    combined$group <- paste0(combined$sample, ":", combined$status)
    combined.wider <- combined %>% select(-sample, -status, -correlation, -ALT) %>% pivot_wider(names_from = cell_type_tissue, values_from = rank) %>% as.data.frame()
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
        # Perform PCA on training data
        train.pca <- prcomp(traindt, scale. = TRUE)  # Scale the data
        
        # Extract top 10 principal components for training
        train_pcs <- train.pca$x[, 1:30]
        
        # Project the test sample onto the PCA space
        project.test <- predict(train.pca, newdata = testdt)
        
        # Ensure test data is in the correct format (1-row matrix with 10 columns)
        test_pcs <- matrix(project.test[1, 1:30], nrow = 1)
        testgroup <- dtmat_group[x]
        
        # Train SVM on the 10 PCs
        svm.model <- svm(train_pcs, traingroup, class.weights = wts, kernel = "linear", cost = cost, scale = FALSE, decision.values = TRUE)
        
        # Predict using SVM
        svm.pred <- predict(svm.model, test_pcs, decision.values = TRUE)
        numeric_true <- ifelse(testgroup == "Healthy", 2, 1)
        alltruegroup <- c(alltruegroup,numeric_true)
        allpredgroup <- c(allpredgroup,attributes(svm.pred)$decision.values)
        numeric_pred <- ifelse(svm.pred == "Healthy", 2, 1)
        allpredval <- c(allpredval,numeric_pred)
        predres <- rbind(predres, data.frame(orgroup[x],attributes(svm.pred)$decision.values))
        #print(data.frame(orgroup[x],svm.pred))
    }
    
                          
    # Assigne control and case levels to 1 and 2
    alltruegroup <- factor(alltruegroup, levels = c(2, 1), labels = c("Control", "Preeclampsia"))
    allpredval <- factor(allpredval, levels = c(2, 1), labels = c("Control", "Preeclampsia"))
    #export decision values and print performance metrics
    #write.table(predres, paste0(output_dir, "/svm_loo_decision_values.txt", sep="\t", row.names=F, col.names=T, quote=F)
    svmperfm <- roc(response=alltruegroup,predictor=allpredgroup,ci=TRUE)
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

    #plot ROC curve
    pdf("/mnt/DATA3/timo/plots/svm_loo_roc_preeclampsia.pdf")
    plot(svmperfm,print.auc=TRUE, main = paste("ROC ", "Preeclampsia"))
    # Add performance metrics to the plot
    text(x = 0.35, y = 0.1, labels = sprintf(
        "Accuracy:      %.2f\nSensitivity:    %.2f\nSpecificity:    %.2f\nPrecision:      %.2f\nF1 Score:      %.2f",
        round(accuracy,2), round(sensitivity,2), round(specificity, 2), round(precision, 2), round(f1_score, 2)),
        pos = 4, cex = 1.0, col = "black")
    dev.off()
}
