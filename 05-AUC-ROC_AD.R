library(lightgbm)
library(dplyr)
library(caret)
library(pROC)

Best_Params <- read.csv("./Results/Best_params.csv")

data1 <- read.csv("./data/Docking_scores/AD_CHIKV_Tetra.csv")
data2 <- read.csv("./data/Docking_scores/AD_DENV_Tetra.csv")
data3 <- read.csv("./data/Docking_scores/AD_WNV_Tetra.csv")
data4 <- read.csv("./data/Docking_scores/AD_ZIKV_Tetra.csv")

Sequence_prop <- read.csv("./data/Properties/Tetra_prop_n.csv")
names(Sequence_prop)[1] <- "Tetra"

pretreat <- function(data){
  names(data) <- c("Tetra", "X")
  data_sequence <- dplyr::left_join(data, Sequence_prop)
  return(data_sequence)
}

data_sequence1 <- pretreat(data1)
data_sequence2 <- pretreat(data2)
data_sequence3 <- pretreat(data3)
data_sequence4 <- pretreat(data4)

prepare.datasets_HO <- function(x, positive_negative = c(0.2, 0.8)){
  x <- x[x$Tetra != "",]
  x <- na.omit(x)
  x[,2] <- as.numeric(x[,2])
  x <- x[order(x[,2]),]
  
  xn <- x[,-c(1:2)]
  rownames(xn) <- x$Tetra
  
  x1 <- replicate(nrow(xn)*positive_negative[1],1)
  x2 <- replicate(nrow(xn)*positive_negative[2],0)
  
  y <- data.frame(y = c(x1, x2))
  y <- rbind(y, data.frame(y = rep(0, nrow(xn) - nrow(y))))
  
  XYn <- as.data.frame(cbind(y, xn))
  
  names(XYn)[1] <- "TARGET"
  return(XYn)
}


train_model <- function(X_train, X_test, y_train, y_test, scale_pos_weight, learning_rate, num_leaves, feature_fraction, pos_bagging_fraction, neg_bagging_fraction, bagging_freq, min_data_in_leaf, max_depth) {
  
  # Set up the parameters
  params <- list(
    objective = "binary",
    metric = "auc",
    boosting = "gbdt",
    scale_pos_weight = scale_pos_weight,
    learning_rate = learning_rate,
    num_leaves = round(num_leaves),
    feature_fraction = feature_fraction,
    pos_bagging_fraction = pos_bagging_fraction,
    neg_bagging_fraction = neg_bagging_fraction,
    bagging_freq = round(bagging_freq),
    min_data_in_leaf = round(min_data_in_leaf),
    max_depth = round(max_depth)
    
  )
  
  dtrain = lgb.Dataset(data = as.matrix(X_train), label = y_train)
  valids = list(validation = lgb.Dataset(data = as.matrix(X_test), label = y_test))
  
  # Train the model
  model <- lgb.train(
    data = dtrain,
    valids = valids,
    params = params,
    nrounds = 5000,
    early_stopping_rounds = 50
  )
  
  return(model)
}

# Train set size
INDEX <- 0.01 

# Positive "best-performing" group
pos_value <- .2 
neg_value <- 1 - pos_value

data_seq1 <- prepare.datasets_HO(data_sequence1, positive_negative = c(pos_value, neg_value))
data_seq2 <- prepare.datasets_HO(data_sequence2, positive_negative = c(pos_value, neg_value))
data_seq3 <- prepare.datasets_HO(data_sequence3, positive_negative = c(pos_value, neg_value))
data_seq4 <- prepare.datasets_HO(data_sequence4, positive_negative = c(pos_value, neg_value))

# Create a list of data sequences
data_sequences <- list(data_seq1, data_seq2, data_seq3, data_seq4)

ROC <- list()

models <- list()
cmDataFrame <- list()
counter <- 1
for (data_seq_idx in 1:length(data_sequences)) {
  
  data_seq <- data_sequences[[data_seq_idx]]
  
  ROC[[data_seq_idx]] <- list()
  
  for (i in 1:3) {
    train_idx <- sample(1:nrow(data_seq), INDEX * nrow(data_seq))
    
    train_data <- data_seq[train_idx, ]
    test_data <- data_seq[-train_idx, ]
    
    X_train <- as.matrix(train_data[, !(names(train_data) %in% c("TARGET", "Tetra"))]) 
    y_train <- train_data$TARGET
    
    X_test <- test_data[, !(names(test_data) %in% c("TARGET", "Tetra"))]
    y_test <- test_data$TARGET
    
    model_results <- train_model(X_train = X_train,
                                 X_test = X_test,
                                 y_train = y_train,
                                 y_test = y_test,
                                 scale_pos_weight = Best_Params$scale_pos_weight, 
                                 learning_rate = Best_Params$learning_rate, 
                                 num_leaves = Best_Params$num_leaves, 
                                 feature_fraction = Best_Params$feature_fraction, 
                                 pos_bagging_fraction = Best_Params$pos_bagging_fraction, 
                                 neg_bagging_fraction = Best_Params$neg_bagging_fraction, 
                                 bagging_freq = Best_Params$bagging_freq, 
                                 min_data_in_leaf = Best_Params$min_data_in_leaf, 
                                 max_depth = Best_Params$max_depth)
    
    # Predict probabilities using the trained model
    y_pred <- predict(model_results, as.matrix(X_test))
    
    ROC[[data_seq_idx]][[i]] <- roc(response = y_test, predictor = y_pred)
    
    counter <- counter + 1
  }
  models[[data_seq_idx]] <- model_results
  cmDataFrame[[data_seq_idx]] <- data.frame(reference = as.factor(y_test), Pred = as.factor(y_pred), Pred_50 = as.factor(round(y_pred)))
}
names(models) <- c("AD_CHIKV", "AD_DENV", "AD_WNV", "AD_ZIKV")
saveRDS(ROC, "./Results/ROCs/ROC-AUC-AD_list.RDs")

saveRDS.lgb.Booster(object = models[[1]], file = "./Results/LightGBM_MODELS/model_AD_CHIKV.RDs")
saveRDS.lgb.Booster(object = models[[2]], file = "./Results/LightGBM_MODELS/model_AD_DENV.RDs")
saveRDS.lgb.Booster(object = models[[3]], file = "./Results/LightGBM_MODELS/model_AD_WNV.RDs")
saveRDS.lgb.Booster(object = models[[4]], file = "./Results/LightGBM_MODELS/model_AD_ZIKV.RDs")

# Plots
ROC <- readRDS("./Results/ROCs/ROC-AUC-AD_list.RDs")
names(ROC) <- c("AD_ZIKV", "AD_ZIKV3", "AD_DENV", "AD_CHIKV", "AD_WNV")
ROC <- ROC[c(4,3,5,1)]

# Plot the ROC curves
mean_auc_values <- numeric(length(ROC))

# Loop through the datasets
for (dataset_idx in 1:length(ROC)) {
  auc_values <- numeric()
  for (i in 1:length(ROC[[dataset_idx]])) {
    roc_obj <- ROC[[dataset_idx]][[i]]
    auc_values[i] <- roc_obj$auc
  }
  # Calculate the mean AUC for this dataset
  mean_auc_values[dataset_idx] <- round(mean(auc_values),2)
}

# Define colors for different datasets
colors <- c("red", "blue", "purple", "orange")
datasets <- c("AD_CHIKV", "AD_DENV", "AD_WNV", "AD_ZIKV")

# Plot ROC curves for each dataset
png(filename = "./Results/Figures/Figure_3.b.png", width = 3500, height = 3200, res = 600)
par(lwd = 1.5, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1, 
    mar = c(5, 4, 4, 2) + 1,   # Adjust the size of margins
    mgp = c(3, 1, 0))  # Move the axis labels away from the plot

for (i in 1:length(ROC)) {
  plot.roc(ROC[[i]][[1]], col = colors[i], add = i != 1, 
           print.thres = T, 
           print.thres.adj = -1000,
           print.thres.pch = 8,
           legacy.axes = T, xlab = "FPR",
           ylab = "TPR")
}

# Add legend
legend("right", legend = paste0(datasets, ": AUC = ", round(mean_auc_values, 2)), 
       col = colors, lty = 1, lwd = 1.5, cex = .7)

thresholds <- data.frame(sapply(ROC, function(x) round(coords(x[[1]], "best"),2)))
thresholds$X1 <- unlist(thresholds$X1)
thresholds$X2 <- unlist(thresholds$X2)
thresholds$X3 <- unlist(thresholds$X3)
thresholds$X4 <- unlist(thresholds$X4)


names(thresholds) <- c("AD_CHIKV", "AD_DENV", "AD_WNV", "AD_ZIKV")

thresholds_transposed <- as.data.frame(t(thresholds))
#write.csv(thresholds, "./Results/ROCs/Best_Thresholds_AD.csv")

colnames(thresholds_transposed) <- c("Threshold", "Specificity", "Sensitivity")

legend("bottomright", legend = paste0("Th = ", thresholds_transposed$Threshold, 
                                      "; TPR = ", thresholds_transposed$Sensitivity, 
                                      "; TNR = ", thresholds_transposed$Specificity),
       col = colors, lty = 1, lwd = 1.5, cex = .7, pch = 8)

dev.off()


ROC_CHIKV <- data.frame(Sensitivity = ROC[[1]][[1]]$sensitivities, Specificity = ROC[[1]][[1]]$specificities, Thresholds = ROC[[1]][[1]]$thresholds)
ROC_DENGV <- data.frame(Sensitivity = ROC[[2]][[1]]$sensitivities, Specificity = ROC[[2]][[1]]$specificities, Thresholds = ROC[[2]][[1]]$thresholds)
ROC_WNV   <- data.frame(Sensitivity = ROC[[3]][[1]]$sensitivities, Specificity = ROC[[3]][[1]]$specificities, Thresholds = ROC[[3]][[1]]$thresholds)
ROC_ZIKV  <- data.frame(Sensitivity = ROC[[4]][[1]]$sensitivities, Specificity = ROC[[4]][[1]]$specificities, Thresholds = ROC[[4]][[1]]$thresholds)

write.csv(ROC_CHIKV, "./Results/ROCs/AD_CHIKV.csv", row.names = F)
write.csv(ROC_DENV, "./Results/ROCs/AD_DENV.csv", row.names = F)
write.csv(ROC_WNV, "./Results/ROCs/AD_WNV.csv", row.names = F)
write.csv(ROC_ZIKV, "./Results/ROCs/AD_ZIKV1.csv", row.names = F)



