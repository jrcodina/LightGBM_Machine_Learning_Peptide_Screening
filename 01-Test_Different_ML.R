#Test different algorithms
library(caret)
library(MLmetrics)
library(gbm)

data <- read.csv("./data/Properties/prop_WNV_tetra_S6n.csv", header= TRUE, sep=",")

prepare.datasets_HO <- function(x, positive_negative = c(0.2, 0.8)){
  x <- x[x$Tetra != "",]
  x <- na.omit(x)
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

data_seq6 <- prepare.datasets_HO(data)

data_seq6$TARGET <- as.factor(data_seq6$TARGET)
levels(data_seq6$TARGET) <- make.names(levels(data_seq6$TARGET)) #Needed for some reason 

ctrl <- trainControl(method = "cv",  classProbs = TRUE, p =0.01, number = 5, savePredictions = "all", returnResamp = "all",  summaryFunction = twoClassSummary)

# SPLS - Sparse Partial Least Square classification
t0 <- Sys.time()
fplsrA <- train(TARGET~., data=data_seq6, method = "spls", trControl = ctrl)
t1 <- Sys.time()
dt_spls <- t1-t0

cm_spls <- confusionMatrix(reference = fplsrA$pred$obs, fplsrA$pred$pred, positive = "X1")
VI_SPLS <- varImp(fplsrA)

# SVM - Support Vector Machine classification
t0 <- Sys.time()
fsvmA <- train(TARGET~., data=data_seq6, method = "svmRadial", trControl = ctrl)
t1 <- Sys.time()
dt_svmRadial <- t1-t0

cm_svm <- confusionMatrix(reference = fsvmA$pred$obs, fsvmA$pred$pred, positive = "X1")
VI_SVM <- varImp(fsvmA)

# RF - Random Forest classification
t0 <- Sys.time()
fRFA <- train(TARGET~., data=data_seq6, method = "rf", trControl = ctrl)
t1 <- Sys.time()
dt_rf <- t1-t0

cm_rf <- confusionMatrix(reference = fRFA$pred$obs, fRFA$pred$pred, positive = "X1")
VI_RF <- varImp(fRFA)

# GBM - Gradient Boost Machine classification
t0 <- Sys.time()
fGBMA <- train(TARGET~., data=data_seq6, method = "gbm", trControl = ctrl)
t1 <- Sys.time()
dt_gbm <- t1-t0

cm_gbm <- confusionMatrix(reference = fGBMA$pred$obs, fGBMA$pred$pred, positive = "X1")
VI_gbm <- varImp(fGBMA)

# RPART - Recursive Partitioning and Regression Tree classification
t0 <- Sys.time()
fRpaA <- train(TARGET~., data=data_seq6, method = "rpart", trControl = ctrl)
t1 <- Sys.time()
dt_rpart <- t1-t0

cm_rpart <- confusionMatrix(reference = fRpaA$pred$obs, fRpaA$pred$pred, positive = "X1")
VI_rpart <- varImp(fRpaA)

# KNN - K-Nearest Neighbors classification
t0 <- Sys.time()
fKNNA <- train(TARGET~., data=data_seq6, method = "knn", trControl = ctrl)
t1 <- Sys.time()
dt_knn <- t1-t0

cm_knn <- confusionMatrix(reference = fKNNA$pred$obs, fKNNA$pred$pred, positive = "X1")
VI_KNN <- varImp(fKNNA)

# Naive_bayes - Naive Bayes classification
t0 <- Sys.time()
fNBA <- train(TARGET~., data=data_seq6, method = "naive_bayes", trControl = ctrl)
t1 <- Sys.time()
dt_naive_bayes <- t1-t0

cm_NB <- confusionMatrix(reference = fNBA$pred$obs, fNBA$pred$pred, positive = "X1")
VI_NB <- varImp(fNBA)

# avNNet - Neural Networks Using Model Averaging
t0 <- Sys.time()
fAvNetA <- train(TARGET~., data=data_seq6, method = "avNNet", trControl = ctrl)
t1 <- Sys.time()
dt_avNNet <- t1-t0

cm_avNNET <- confusionMatrix(reference = fAvNetA$pred$obs, fAvNetA$pred$pred, positive = "X1")
VI_avNNET <- varImp(fAvNetA)

# NNET - Neural Network
t0 <- Sys.time()
fNNetA <- train(TARGET~., data=data_seq6, method = "nnet", trControl = ctrl)
t1 <- Sys.time()
dt_nnet <- t1-t0

cm_nnet <- confusionMatrix(reference = fNNetA$pred$obs, fNNetA$pred$pred, positive = "X1")
VI_NNET <- varImp(fNNetA)


library(lightgbm)
data_seq6 <- prepare.datasets_HO(data)

train_idx <- sample(1:nrow(data_seq6), 0.01 * nrow(data_seq6))

train_data <- data_seq6[train_idx, ]
test_data <- data_seq6[-train_idx, ]

X_train <- as.matrix(train_data[, !(names(train_data) %in% c("TARGET", "Tetra"))]) 
y_train <- train_data$TARGET

X_test <- test_data[, !(names(test_data) %in% c("TARGET", "Tetra"))]
y_test <- test_data$TARGET


# Parameters for LightGBM, similar to caret's trainControl
params <- list(objective = "binary", 
               metric = "auc")

# Time and run LightGBM with 5-fold cross-validation
t0 <- Sys.time()
cv_lgb <- lgb.cv(data = as.matrix(X_train),
                 label = y_train,
                 params = params,
                 nrounds = 100, 
                 nfold = 5,
                 early_stopping_rounds = 20)

# Train a full model with the best iteration
bst <- lgb.train(params = params,
                 data = as.matrix(X_train),
                 label = y_train,
                 nrounds = best_iter)

# Predict and compute confusion matrix
preds <- predict(bst, as.matrix(X_test))
preds <- round(preds)

t1 <- Sys.time()
dt_lgbm <- t1-t0

cm_lgbm <- confusionMatrix(reference = as.factor(y_test), as.factor(preds), positive = "1")

# Feature importance
VI_lgbm <- lgb.importance(bst)


# Compare all

all_algo <- readRDS("../algorithm_comparison.RDS")

time <- data.frame(Method = c("LGBM", "Naive_Bayes", "Rpart", "GBM", "NNET", "avNNET", "KNN", "RF", "SVM-Radial"),
                   Type = c("Tree-Based", "Linear-Based", "Tree-Based", "Tree-Based", "Neural-Network-Based", "Neural-Network-Based", "Kernel-Based", "Tree-Based", "Kernel-Based"),
                   Time = c(dt_lgbm, dt_naive_bayes, dt_rpart, dt_gbm, dt_nnet, dt_avNNet, dt_knn, dt_rf, dt_svmRadial), 
                   Accuracy = c(cm_lgbm[["overall"]][["Accuracy"]], cm_NB[["overall"]][["Accuracy"]], 
                                cm_rpart[["overall"]][["Accuracy"]], cm_gbm[["overall"]][["Accuracy"]],
                                cm_nnet[["overall"]][["Accuracy"]], cm_avNNET[["overall"]][["Accuracy"]],
                                cm_knn[["overall"]][["Accuracy"]], cm_rf[["overall"]][["Accuracy"]],
                                cm_svm[["overall"]][["Accuracy"]]))

all_algorithms <- list(dataset = data_seq6, 
                       control = ctrl,
                       LGBM = list(confussionMatrix = cm_lgbm, model = bst, Variant.importance = VI_lgbm, time = dt_lgbm),
                       Naive_Bayes = list(confussionMatrix = cm_NB, model = fNBA, Variant.importance = VI_NB, time = dt_naive_bayes),
                       Rpart = list(confussionMatrix = cm_rpart, model = fRpaA, Variant.importance = VI_rpart, time = dt_rpart),
                       GBM = list(confussionMatrix = cm_gbm, model = fGBMA, Variant.importance = VI_gbm, time = dt_gbm),
                       NNET = list(confussionMatrix = cm_nnet, model = fNNetA, Variant.importance = VI_NNET, time = dt_nnet),
                       avNNET = list(confussionMatrix = cm_avNNET, model = fAvNetA, Variant.importance = VI_avNNET, time = dt_avNNet),
                       KNN = list(confussionMatrix = cm_knn, model = fKNNA, Variant.importance = VI_KNN, time = dt_knn),
                       RF = list(confussionMatrix = cm_rf, model = fRFA, Variant.importance = VI_RF, time = dt_rf),
                       SVM_Radial = list(confussionMatrix = cm_svm, model = fsvmA, Variant.importance = VI_SVM, time = dt_svmRadial))
saveRDS(all_algorithms, "./Supplementary/algorithm_comparison.RDS")


results_df <- data.frame(Names = character(), F1 = numeric(), Time = numeric())

for (i in 3:11) {
  name <- names(all_alg)[i]
  F1 <- round(all_alg[[i]][["confussionMatrix"]][["byClass"]][["F1"]],2)
  time <- as.numeric(all_alg[[i]][["time"]], units = "mins")
  
  # Store results in data.frame
  results_df[i - 2,] <- c(name, F1, time)
}

write.csv(results_df, "./Results/Table_1.csv", row.names = F)
