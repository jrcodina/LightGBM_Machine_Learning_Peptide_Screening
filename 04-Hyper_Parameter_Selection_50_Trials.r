#dt = 2.52 h
t0 <- Sys.time()

library(lightgbm)
library(dplyr)
library(rBayesianOptimization)
library(data.table)
library(foreach)
library(doParallel)

data <- read.csv("./data/Docking_scores/WNV_Tetra.csv")
names(data) <- c("Tetra", "X")

Sequence_prop <- read.csv("./data/Properties/Tetra_prop_n.csv")
names(Sequence_prop)[1] <- "Tetra"

data_sequence <- dplyr::left_join(data, Sequence_prop)

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


train_model <- function(scale_pos_weight, learning_rate, num_leaves, feature_fraction, pos_bagging_fraction, neg_bagging_fraction, bagging_freq, min_data_in_leaf, max_depth) {
  
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
  
  # Train the model
  model <- lgb.train(
    data = dtrain,
    valids = list(validation = lgb.Dataset(data = as.matrix(X_test), label = y_test)),
    params = params,
    nrounds = 5000,
    early_stopping_rounds = 50
  )
  
  first_pred <- predict(model, as.matrix(X_test))
  
  # Evaluate the model
  all_prop <- data_seq[,-1]
  pred <- predict(model, as.matrix(all_prop), reshape = T)
  
  pred_y = round(first_pred)
  
  if (all(pred_y == 0)) {
    pred_y[1] <- 1
  }
  
  real_vs_pred <- data.frame(Real = y_test, Pred = pred_y)
  
  cm <- caret::confusionMatrix(reference = as.factor(real_vs_pred$Real), as.factor(real_vs_pred$Pred), positive = "1")
  
  F1 <- cm[["byClass"]][["F1"]]
  F1 <- ifelse(is.na(F1), 0, F1)
  
  return(list(Score = F1, Pred = first_pred))
}

# Positive "best-performing" group
pos_value <- .2
neg_value <- 1 - pos_value

# Train set size
INDEX <- 0.01

data_seq <- prepare.datasets_HO(data_sequence, c(pos_value, neg_value))

# Initialize a list to store the results
results <- list()

# Register the parallel backend
cl <- makeCluster(8)
registerDoParallel(cl)

# Loop through the positive values
results <- foreach(i = 1:50, .combine = "list", .packages = c("lightgbm", "caret", "rBayesianOptimization", "dplyr"), .export = c("prepare.datasets_HO", "data_sequence", "INDEX", "train_model")) %dopar% {
  
  train_idx <- sample(1:nrow(data_seq), INDEX * nrow(data_seq))
  
  train_data <- data_seq[train_idx, ]
  test_data <- data_seq[-train_idx, ]
  
  X_train <- as.matrix(train_data[, !(names(train_data) %in% c("TARGET", "Tetra"))]) 
  y_train <- train_data$TARGET
  
  X_test <- test_data[, !(names(test_data) %in% c("TARGET", "Tetra"))]
  y_test <- test_data$TARGET
  
  bounds <- list(
    scale_pos_weight = c(1L, 50L),
    learning_rate = c(0.001, 0.9),
    num_leaves = c(8L, 31L),
    feature_fraction = c(0.1, 1),
    pos_bagging_fraction = c(0.1, 0.99),
    neg_bagging_fraction = c(0.1, 0.99),
    bagging_freq = c(1L, 30L),
    min_data_in_leaf = c(5L, 90L),
    max_depth = c(1L, 10L)
  )
  
  opt_result_Sequence2 <- BayesianOptimization(
    FUN = train_model,
    bounds = bounds,
    init_points = 25,
    n_iter = 15,
    acq = "ucb",
    kappa = 0.5,
    verbose = TRUE
  )
  
  results[[i]] <- list(Positives = pos_value, INDEX = INDEX, F1 = opt_result_Sequence2$Best_Value, Params = opt_result_Sequence2$Best_Par)
}


stopCluster(cl)

b <- unlist(results, recursive = T)
num_cols <- 12

c <- data.frame(Positive =             b[seq(1, length(b-(num_cols+1)), num_cols)],
                INDEX =                b[seq(2, length(b-(num_cols+2)), num_cols)],
                F1 =                   b[seq(3, length(b-(num_cols+3)), num_cols)], 
                scale_pos_weight =     b[seq(4, length(b-(num_cols+4)), num_cols)], 
                learning_rate =        b[seq(5, length(b-(num_cols+5)), num_cols)], 
                num_leaves =           b[seq(6, length(b-(num_cols+6)), num_cols)],
                feature_fraction =     b[seq(7, length(b-(num_cols+7)), num_cols)],
                pos_bagging_fraction = b[seq(8, length(b-(num_cols+8)), num_cols)],
                neg_bagging_fraction = b[seq(9, length(b-(num_cols+9)), num_cols)],
                bagging_freq =         b[seq(10, length(b-(num_cols+10)), num_cols)],
                min_data_in_leaf =     b[seq(11, length(b-(num_cols+11)), num_cols)],
                max_depth =            b[seq(12, length(b-(num_cols+12)), num_cols)])


t1 <- Sys.time()
dt <- t1-t0
dt

best_params <- c[c$F1 == max(c$F1),]

write.csv(c, "./Not_in_paper/Hyperparam_selection.csv", row.names = F)
write.csv(best_params, "./Results/Best_params.csv", row.names = F)
