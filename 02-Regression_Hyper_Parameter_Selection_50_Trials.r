#dt = 2.52 h
t0 <- Sys.time()

library(lightgbm)
library(dplyr)
library(rBayesianOptimization)
library(data.table)
library(foreach)
library(doParallel)

data <- read.csv("./data/Properties/prop_WNV_tetra_S6n.csv")
names(data)[2] <- "TARGET"
data_seq <- data

train_model <- function(learning_rate, num_leaves, feature_fraction, bagging_freq, min_data_in_leaf, max_depth) {
  
  # Set up the parameters
  params <- list(
    objective = "regression",
    metric = "mae",
    boosting = "gbdt",
    learning_rate = learning_rate,
    num_leaves = round(num_leaves),
    feature_fraction = feature_fraction,
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
  
  return(list(Score = model$best_score, Pred = first_pred))
}

# Train set size
INDEX <- 0.01

# Initialize a list to store the results
results <- list()

# Register the parallel backend
cl <- makeCluster(2)
registerDoParallel(cl)

# Loop through the positive values
results <- foreach(i = 1:5, .combine = "list", .packages = c("lightgbm", "caret", "rBayesianOptimization", "dplyr"), .export = c("INDEX", "train_model")) %dopar% {
  
  train_idx <- sample(1:nrow(data_seq), INDEX * nrow(data_seq))
  
  train_data <- data_seq[train_idx, ]
  test_data <- data_seq[-train_idx, ]
  
  X_train <- as.matrix(train_data[, !(names(train_data) %in% c("TARGET", "Tetra"))]) 
  y_train <- train_data$TARGET
  
  X_test <- test_data[, !(names(test_data) %in% c("TARGET", "Tetra"))]
  y_test <- test_data$TARGET
  
  bounds <- list(
    learning_rate = c(0.001, 0.9),
    num_leaves = c(8L, 31L),
    feature_fraction = c(0.1, 1),
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
  
  results[[i]] <- list(MAE = opt_result_Sequence2$Best_Value, Params = opt_result_Sequence2$Best_Par)
}

stopCluster(cl)

b <- unlist(results, recursive = T)

num_cols <- 7

c <- data.frame(MAE =                  b[seq(1, (length(b)-(num_cols+1)), num_cols)], 
                learning_rate =        b[seq(2, length(b-(num_cols+2)), num_cols)], 
                num_leaves =           b[seq(3, length(b-(num_cols+3)), num_cols)],
                feature_fraction =     b[seq(4, length(b-(num_cols+4)), num_cols)],
                bagging_freq =         b[seq(5, length(b-(num_cols+5)), num_cols)],
                min_data_in_leaf =     b[seq(6, length(b-(num_cols+6)), num_cols)],
                max_depth =            b[seq(7, length(b-(num_cols+7)), num_cols)])

best_params <- c[c$MAE == min(c$MAE),]

write.csv(c, "./Not_in_paper/Regression.csv", row.names = F)
write.csv(best_params, "./Results/Best_params_regress.csv", row.names = F)

# Plot

best_params <- read.csv("./Results/Best_params_regress.csv")

train_idx <- sample(1:nrow(data_seq), 0.01 * nrow(data_seq))

train_data <- data_seq[train_idx, ]
test_data <- data_seq[-train_idx, ]

X_train <- as.matrix(train_data[, !(names(train_data) %in% c("TARGET", "Tetra"))]) 
y_train <- train_data$TARGET

X_test <- test_data[, !(names(test_data) %in% c("TARGET", "Tetra"))]
y_test <- test_data$TARGET

params <- list(
  objective = "regression",
  metric = "mae",
  boosting = "gbdt",
  learning_rate = best_params$learning_rate,
  num_leaves = round(best_params$num_leaves),
  feature_fraction = best_params$feature_fraction,
  bagging_freq = round(best_params$bagging_freq),
  min_data_in_leaf = round(best_params$min_data_in_leaf),
  max_depth = round(best_params$max_depth)
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

real_vs_pred <- data.frame(Tetra = test_data$Tetra, Real = test_data$TARGET, Pred = first_pred)

prepare.datasets <- function(x, positive_negative = c(0.2, 0.8), score.col = 2, prefix){
  x <- x[x$Tetra != "",]
  x <- na.omit(x)
  x <- x[order(x[,score.col]),]
  
  x1 <- replicate(nrow(x)*positive_negative[1],1)
  x2 <- replicate(nrow(x)*positive_negative[2],0)
  
  y <- data.frame(y = c(x1, x2))
  y <- rbind(y, data.frame(y = rep(0, nrow(x) - nrow(y))))
  
  XYn <- as.data.frame(cbind(x, y))
  
  names(XYn)[ncol(x)+1] <- paste0(prefix, ".Binary")
  return(XYn)
}

real_vs_pred <- prepare.datasets(positive_negative = c(0.2, 0.8), real_vs_pred, score.col = 2, prefix="Real")
real_vs_pred <- prepare.datasets(positive_negative = c(0.2, 0.8), real_vs_pred, score.col = 3, prefix="Pred")

cmRegression <- caret::confusionMatrix(reference = as.factor(real_vs_pred$Real.Binary), as.factor(real_vs_pred$Pred.Binary), positive = "1")
df <- real_vs_pred

# Filter data for each subset
lower_bound_limit <- quantile(df$Real, 0.05)
upper_bound_limit <- quantile(df$Real, 0.95)

lower_bound <- real_vs_pred %>% filter(Real < lower_bound_limit)
upper_bound <- real_vs_pred %>% filter(Real > upper_bound_limit)
middle_bound <- real_vs_pred %>% filter(Real >= lower_bound_limit & Real <= upper_bound_limit)

# Calculate RMSE for each subset
rmse <- sqrt(mean((real_vs_pred$Real - real_vs_pred$Pred)^2))

rmse_lower_bound <- sqrt(mean((lower_bound$Real - lower_bound$Pred)^2))
rmse_upper_bound <- sqrt(mean((upper_bound$Real - upper_bound$Pred)^2))
rmse_middle_bound <- sqrt(mean((middle_bound$Real - middle_bound$Pred)^2))

library(ggplot2)
library(ggExtra)

p <- ggplot(real_vs_pred, aes(x = Real, y = Pred)) +
  geom_point(color = "steelblue", alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = lower_bound_limit, linetype = "dotted", color = "gray30", size = 1.3) +
  geom_vline(xintercept = upper_bound_limit, linetype = "dotted", color = "gray30", size = 1.3) +
  annotate(size = 8, "text", x = -7, y = 9.5, label = paste0("RMSE: ", round(rmse_lower_bound, 2))) +
  annotate(size = 8,"text", x = 6.5, y = 9.5, label = paste0("RMSE: ", round(rmse_upper_bound, 2))) +
  annotate(size = 8,"text", x = -0.5, y = 9.5, label = paste0("RMSE: ", round(rmse_middle_bound, 2))) +
  theme_bw() +
  theme(text = element_text(size = 36, color = "black"),
        axis.title = element_text(size = 36, color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=3), 
        axis.text = element_text(size = 36, color = "black"), 
        plot.caption = element_text(size = 20),
        panel.grid = element_line(color = "gray")) +
  labs(x = "Real Score",
       y = "Predicted Score",
       caption = paste0("RMSE (All Segments): ", round(rmse,2))) +
  scale_x_continuous(labels = scales::comma, limits = c(-10, 10)) +
  scale_y_continuous(labels = scales::comma, limits = c(-10, 10))

pm <- ggMarginal(p, type = "histogram")

ggsave(plot = pm, device = "png", paste0("./Supplementary/Figure_S3.png"),width = 10, height = 10, dpi = 1200)


#----------------------------#

# Other plot options

df <- real_vs_pred

ggplot(df, aes(x = Real, y = Pred)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  geom_abline(slope = 1, intercept = 0, color = 'red', linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, hjust = 0.5)) +
  labs(x = "Real Docking Scores",
       y = "Predicted Docking Scores",
       title = "Real vs Predicted Docking Scores",
       fill = "Density")

library(hexbin)

ggplot(df, aes(x = Real, y = Pred)) +
  geom_hex(bins = 50) +  # you can adjust the number of bins as needed
  geom_abline(slope = 1, intercept = 0, color = 'white', linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, hjust = 0.5)) +
  labs(x = "Real Docking Scores",
       y = "Predicted Docking Scores",
       title = "Real vs Predicted Docking Scores",
       fill = "Count")

library(ggExtra)

df <- real_vs_pred

p <- ggplot(df, aes(x = Real, y = Pred)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = 'red', linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, hjust = 0.5)) +
  labs(x = "Real Docking Scores",
       y = "Predicted Docking Scores",
       title = "Real vs Predicted Docking Scores")


ggplot(df, aes(x = Real, y = Pred)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = 'red', linetype = "dashed") +
  stat_density_2d(aes(color = ..level..), geom = "contour", bins = 10) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, hjust = 0.5)) +
  labs(x = "Real Docking Scores",
       y = "Predicted Docking Scores",
       title = "Real vs Predicted Docking Scores",
       color = "Density Level")

ggplot(real_vs_pred, aes(y = Real, x = 1:nrow(real_vs_pred))) +
  geom_line(size = 1, aes(y = Pred), color = "red")+
  geom_line(size = 2) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, hjust = 0.5)) +
  labs(x = "Numer of Peptides",
       y = "Docking Scores",
       color = "Receptor") +
  scale_color_manual(values = colors)+
  scale_x_continuous(labels = scales::comma)

library(ggplot2)

ggplot(real_vs_pred, aes(x = Real, y = Pred)) +
  geom_point(color = "blue") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, hjust = 0.5)) +
  labs(x = "Real Scores",
       y = "Predicted Scores") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma)

library(ggplot2)


ggplot(real_vs_pred, aes(y = Real, x = Pred)) +
  geom_point(color = "steelblue", alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, hjust = 0.5)) +
  labs(x = "Predicted Score",
       y = "Real Score",
       title = "Comparison of Real Scores vs Predicted Scores",
       caption = paste0("Data source: real_vs_pred\nRMSE: ", round(rmse, 2))) +
  scale_x_continuous(labels = scales::comma, limits = c(min(real_vs_pred$Real)-0.5, max(real_vs_pred$Real)+0.5)) +
  scale_y_continuous(labels = scales::comma, limits = c(min(real_vs_pred$Real)-0.5, max(real_vs_pred$Real)+0.5))

