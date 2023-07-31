#dt = 18.21 mins
t0 <- Sys.time()

library(lightgbm)
library(dplyr)
library(tidyr)
library(data.table)


Best_Params <- read.csv("./Results/Best_params.csv")

data1 <- read.csv("./data/Docking_scores/CHIKV_Tetra.csv")
data2 <- read.csv("./data/Docking_scores/DENV_Tetra.csv")
data3 <- read.csv("./data/Docking_scores/WNV_Tetra.csv")
data4 <- read.csv("./data/Docking_scores/ZIKV_Tetra.csv")

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

data_seq1 <- prepare.datasets_HO(data_sequence1)
data_seq2 <- prepare.datasets_HO(data_sequence2)
data_seq3 <- prepare.datasets_HO(data_sequence3)
data_seq4 <- prepare.datasets_HO(data_sequence4)

train_model <- function(reference_data, opt_th, scale_pos_weight, learning_rate, num_leaves, feature_fraction, pos_bagging_fraction, neg_bagging_fraction, bagging_freq, min_data_in_leaf, max_depth) {
  
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
  
  first_pred <- predict(model, as.matrix(X_test))
  
  # Evaluate the model
  all_prop <- reference_data[,-1]
  pred <- predict(model, as.matrix(all_prop), reshape = T)
  
  pred_y = ifelse(pred >= opt_th, 1, 0)
  
  if (all(pred_y == 0)) {
    pred_y[1] <- 1
  }
  
  real_vs_pred <- data.frame(Real = reference_data[,1], Pred = pred_y)
  
  cm <- caret::confusionMatrix(reference = as.factor(real_vs_pred$Real), as.factor(real_vs_pred$Pred), positive = "1")
  
  Precision <- cm[["byClass"]][["Precision"]]
  specificity <- cm[["byClass"]][["Specificity"]]
  sensitivity <- cm[["byClass"]][["Sensitivity"]]
  Accuracy <- cm[["overall"]][["Accuracy"]]
  F1 <- cm[["byClass"]][["F1"]]
  kappa <- cm[["overall"]][["Kappa"]]
  TP <- cm[["table"]]["1","1"]
  FP <- cm[["table"]]["1","0"]
  FN <- cm[["table"]]["0","1"]
  TN <- cm[["table"]]["0","0"]
  
  return(data.frame(Kappa = kappa, Precision = Precision, Specificity = specificity, Sensitivity = sensitivity, TP = TP, FP = FP, FN = FN, TN = TN, F1 = F1, Accuracy = Accuracy))
}

INDEX <- 0.01

pos_value <- .2
neg_value <- 1 - pos_value


# Initialize a list to store the results
results <- data.frame(Kappa = numeric(), Precision = numeric(), Specificity = numeric(), Sensitivity = numeric(), TP = numeric(), FP = numeric(), FN = numeric(), TN = numeric(), F1 = numeric(), Accuracy = numeric())

# Create a list of data sequences
data_sequences <- list(data_seq1, data_seq2, data_seq3, data_seq4, data_seq5, data_seq6)

# Loop through each data sequence
for (j in 1:length(data_sequences)) {
  data_seq<- data_sequences[[j]]
  opt_th <- optimal_th[1,j]
  # Test the data sequence with 100 iterations
  for (i in 1:100) {
    train_idx <- sample(1:nrow(data_seq), INDEX * nrow(data_seq))
    
    train_data <- data_seq[train_idx, ]
    test_data <- data_seq[-train_idx, ]
    
    X_train <- as.matrix(train_data[, !(names(train_data) %in% c("TARGET", "Tetra"))]) 
    y_train <- train_data$TARGET
    
    X_test <- test_data[, !(names(test_data) %in% c("TARGET", "Tetra"))]
    y_test <- test_data$TARGET
    
    metrics <- train_model(reference_data = data_seq,
                                 opt_th = opt_th,
                                 scale_pos_weight = Best_Params$scale_pos_weight, 
                                 learning_rate = Best_Params$learning_rate, 
                                 num_leaves = Best_Params$num_leaves, 
                                 feature_fraction = Best_Params$feature_fraction, 
                                 pos_bagging_fraction = Best_Params$pos_bagging_fraction, 
                                 neg_bagging_fraction = Best_Params$neg_bagging_fraction, 
                                 bagging_freq = Best_Params$bagging_freq, 
                                 min_data_in_leaf = Best_Params$min_data_in_leaf, 
                                 max_depth = Best_Params$max_depth)
    
    results <- rbind(results, metrics)
  }
}


t1 <- Sys.time()
dt <- t1-t0
dt

results$DataSet <- c(rep("CHIKV", 100), rep("DENV", 100), rep("WNV", 100), rep("ZIKV", 100))

df <- results

df_Accuracy <- df %>%
  group_by(DataSet) %>%
  summarize(Mean = mean(Accuracy, na.rm = TRUE),
            SD = sd(Accuracy, na.rm = TRUE))

df_Specificity <- df %>%
  group_by(DataSet) %>%
  summarize(Mean = mean(Specificity, na.rm = TRUE),
            SD = sd(Specificity, na.rm = TRUE))

df_Sensitivity <- df %>%
  group_by(DataSet) %>%
  summarize(Mean = mean(Sensitivity, na.rm = TRUE),
            SD = sd(Sensitivity, na.rm = TRUE))

df_F1 <- df %>%
  group_by(DataSet) %>%
  summarize(Mean = mean(F1, na.rm = TRUE),
            SD = sd(F1, na.rm = TRUE))

df_Kappa <- df %>%
  group_by(DataSet) %>%
  summarize(Mean = mean(Kappa, na.rm = TRUE),
            SD = sd(Kappa, na.rm = TRUE))

df_Precision <- df %>%
  group_by(DataSet) %>%
  summarize(Mean = mean(Precision, na.rm = TRUE),
            SD = sd(Precision, na.rm = TRUE))

a <- dplyr::left_join(suffix = c("Acc", "Speci"),df_Accuracy, df_Specificity, by = "DataSet")
b <- dplyr::left_join(suffix = c("Sensi", "F1"),df_Sensitivity, df_F1, by = "DataSet")
c <- dplyr::left_join(suffix = c("Kappa", "PPV"),df_Kappa, df_Precision, by = "DataSet")

Table_3.1 <- cbind(a,b,c)
Table_3.1 <- Table_3.1[,c(1, 2, 3, 7, 8, 4,5, 9, 10)]
Table_3.1 <- Table_3.1[c(4,5,6,2,1,3),]
Table_3.1[,-1] <- apply(Table_3.1[,-1], 2, function(x) round(x, 3))

write.csv(Table_3.1, "./Results/Table_3.1.csv", row.names = F)
write.csv(df, "./Results/Table_3.1_Full_Data.csv", row.names = F)