library(lightgbm)
library(dplyr)

Best_Params <- read.csv("./Results/Best_params.csv")

data1 <- read.csv("./data/Docking_scores/CHIKV_Tetra.csv")
data2 <- read.csv("./data/Docking_scores/DENV_Tetra.csv")
data3 <- read.csv("./data/Docking_scores/WNV_Tetra.csv")
data4 <- read.csv("./data/Docking_scores/ZIKV_Tetra.csv")

data1AD <- read.csv("./data/Docking_scores/AD_CHIKV_Tetra.csv")
data2AD <- read.csv("./data/Docking_scores/AD_DENV_Tetra.csv")
data3AD <- read.csv("./data/Docking_scores/AD_WNV_Tetra.csv")
data4AD <- read.csv("./data/Docking_scores/AD_ZIKV_Tetra.csv")

Sequence_prop <- read.csv("./data/Properties/Tetra_prop_n.csv")
names(Sequence_prop)[1] <- "Tetra"

pretreat <- function(data){
  data <- data[,1:2]
  names(data) <- c("Tetra", "X")
  data <- data[order(data$X),]
  data_sequence <- dplyr::left_join(data, Sequence_prop)
  return(data_sequence)
}

data_sequence1 <- pretreat(data1)
data_sequence2 <- pretreat(data2)
data_sequence3 <- pretreat(data3)
data_sequence4 <- pretreat(data4)

data_sequence1AD <- pretreat(data1AD)
data_sequence2AD <- pretreat(data2AD)
data_sequence3AD <- pretreat(data3AD)
data_sequence4AD <- pretreat(data4AD)

prepare.datasets_HO <- function(x, positive_negative = c(0.2, 0.8)){
  x <- x[x$Tetra != "",]
  x <- na.omit(x)
  x <- x[order(x$X),]
  
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

data_seq1AD <- prepare.datasets_HO(data_sequence1AD)
data_seq2AD <- prepare.datasets_HO(data_sequence2AD)
data_seq3AD <- prepare.datasets_HO(data_sequence3AD)
data_seq4AD <- prepare.datasets_HO(data_sequence4AD)

data_sequences <- list(data_seq1,   data_seq2,   data_seq3,   data_seq4,
                       data_seq1AD, data_seq2AD, data_seq3AD, data_seq4AD)

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

models <- list()

for (i in 1:length(data_sequences)) {
  data_seq <- data_sequences[[i]]
  
  train_idx <- sample(1:nrow(data_seq), 0.01 * nrow(data_seq))
  
  train_data <- data_seq[train_idx, ]
  test_data <- data_seq[-train_idx, ]
  
  X_train <- as.matrix(train_data[, !(names(train_data) %in% c("TARGET", "Tetra"))]) 
  y_train <- train_data$TARGET
  
  X_test <- test_data[, !(names(test_data) %in% c("TARGET", "Tetra"))]
  y_test <- test_data$TARGET
  
  models[[i]] <- train_model(X_train = X_train,
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
}

names(models) <- c("CHIKV",    "DENV",    "WNV",    "ZIKV",
                   "CHIKV_AD", "DENV_AD", "WNV_AD", "ZIKV_AD")

pred_data1 <- data.frame(Tetra = rownames(data_seq1), Pred = predict(models$CHIKV, as.matrix(data_seq1[,-1])), Real = data_seq1$TARGET, score = na.omit(data_sequence1$X))
pred_data2 <- data.frame(Tetra = rownames(data_seq2), Pred = predict(models$DENV, as.matrix(data_seq2[,-1])), Real = data_seq2$TARGET, score = na.omit(data_sequence2$X))
pred_data3 <- data.frame(Tetra = rownames(data_seq3), Pred = predict(models$WNV, as.matrix(data_seq3[,-1])), Real = data_seq3$TARGET, score = na.omit(data_sequence3$X))
pred_data4 <- data.frame(Tetra = rownames(data_seq4), Pred = predict(models$ZIKV, as.matrix(data_seq4[,-1])), Real = data_seq4$TARGET, score = na.omit(data_sequence4$X))

pred_data1AD <- data.frame(Tetra = rownames(data_seq1AD), Pred = predict(models$CHIKV_AD, as.matrix(data_seq1AD[,-1])), Real = data_seq1AD$TARGET, score = na.omit(data_sequence1AD$X))
pred_data2AD <- data.frame(Tetra = rownames(data_seq2AD), Pred = predict(models$DENV_AD, as.matrix(data_seq2AD[,-1])), Real = data_seq2AD$TARGET, score = na.omit(data_sequence2AD$X))
pred_data3AD <- data.frame(Tetra = rownames(data_seq3AD), Pred = predict(models$WNV_AD, as.matrix(data_seq3AD[,-1])), Real = data_seq3AD$TARGET, score = na.omit(data_sequence3AD$X))
pred_data4AD <- data.frame(Tetra = rownames(data_seq4AD), Pred = predict(models$ZIKV_AD, as.matrix(data_seq4AD[,-1])), Real = data_seq4AD$TARGET, score = na.omit(data_sequence4AD$X))

# STRUCTURAL ANALYSIS

# Count the occurrences of each amino acid in the current position
count_aa <- function(Title, n = 4) {
  aa_counts <- data.frame(matrix(nrow = n, ncol = 20))
  colnames(aa_counts) <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  
  for (w in 1:n) {
    position_aa <- substr(Title, w, w)
    aa_counts[w,] <- table(factor(position_aa, levels = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")))
  }
  return(aa_counts)
}


# Then we transform to -1, 0, 1 matrix based on a probability of appearance of each peptide in each position. The probability threshold is set in the argument 'threshold'.
transform_aa <- function(x, threshold = 1.65, group = "Best", pred.col, percentage = 0.05) {
  
  total <- count_aa(x$Tetra)
  
  total_sample <- nrow(x)
  
  if (group == "Best") {
    sample_size = length(x$Tetra[1:(nrow(x)*percentage)])
    selected = count_aa(x$Tetra[1:(nrow(x)*percentage)])
  } 
  
  if (group == "Worst") {
    sample_size = length(x$Tetra[(nrow(x)-nrow(x)*percentage):nrow(x)])
    selected = count_aa (x$Tetra[(nrow(x)-nrow(x)*percentage):nrow(x)])
  }
  
  probability_matrix <- ((selected/sample_size) - ((total-selected) / (total_sample-sample_size))) / (sqrt((total/total_sample)*(1-total/total_sample))* sqrt(1/sample_size + 1/(total_sample-sample_size))) # To get the probability
  
  transformed <- ifelse(probability_matrix > threshold, 1, ifelse(probability_matrix < -threshold, -1, 0)) # TO transform the probability in -1, 0, 1 based on threshold
  
  return(transformed)
}

transformed.best.1 <- transform_aa(pred_data1, threshold = 1.65, group = "Best", percentage = 0.05)
transformed.best.2 <- transform_aa(pred_data2, threshold = 1.65, group = "Best", percentage = 0.05)
transformed.best.3 <- transform_aa(pred_data3, threshold = 1.65, group = "Best", percentage = 0.05)
transformed.best.4 <- transform_aa(pred_data4, threshold = 1.65, group = "Best", percentage = 0.05)

transformed.best.1AD <- transform_aa(pred_data1AD, threshold = 1.65, group = "Best", percentage = 0.05)
transformed.best.2AD <- transform_aa(pred_data2AD, threshold = 1.65, group = "Best", percentage = 0.05)
transformed.best.3AD <- transform_aa(pred_data3AD, threshold = 1.65, group = "Best", percentage = 0.05)
transformed.best.4AD <- transform_aa(pred_data4AD, threshold = 1.65, group = "Best", percentage = 0.05)

transformed.worst.1 <- transform_aa(pred_data1, threshold = 1.65,group = "Worst", percentage = 0.05)
transformed.worst.2 <- transform_aa(pred_data2, threshold = 1.65,group = "Worst", percentage = 0.05)
transformed.worst.3 <- transform_aa(pred_data3, threshold = 1.65,group = "Worst", percentage = 0.05)
transformed.worst.4 <- transform_aa(pred_data4, threshold = 1.65,group = "Worst", percentage = 0.05)

transformed.worst.1AD <- transform_aa(pred_data1AD, threshold = 1.65,group = "Worst", percentage = 0.05)
transformed.worst.2AD <- transform_aa(pred_data2AD, threshold = 1.65,group = "Worst", percentage = 0.05)
transformed.worst.3AD <- transform_aa(pred_data3AD, threshold = 1.65,group = "Worst", percentage = 0.05)
transformed.worst.4AD <- transform_aa(pred_data4AD, threshold = 1.65,group = "Worst", percentage = 0.05)

# Select the ones that are 1 for best and 0 or -1 in worst and expand the grid to all possible combinations
best.vs.worst <- function(transformed.best, transformed.worst){
  x <- as.data.frame(t(transformed.best))
  y <- as.data.frame(t(transformed.worst))
  
  poss <- c()
  for (j in 1:length(x)) {
    poss[[j]] <-  dplyr::setdiff(rownames(dplyr::filter(x, x[,j] == 1)),
                                 rownames(dplyr::filter(y, y[,j] == 1)))
  }
  
  best.pred <- expand.grid(poss)
  best.pred <- data.frame(Title = apply(best.pred, 1, paste, collapse = ""))
  return(best.pred)
}

SA_data1 <- best.vs.worst(transformed.best.1, transformed.worst.1)
SA_data2 <- best.vs.worst(transformed.best.2, transformed.worst.2)
SA_data3 <- best.vs.worst(transformed.best.3, transformed.worst.3)
SA_data4 <- best.vs.worst(transformed.best.4, transformed.worst.4)

SA_data1AD <- best.vs.worst(transformed.best.1AD, transformed.worst.1AD)
SA_data2AD <- best.vs.worst(transformed.best.2AD, transformed.worst.2AD)
SA_data3AD <- best.vs.worst(transformed.best.3AD, transformed.worst.3AD)
SA_data4AD <- best.vs.worst(transformed.best.4AD, transformed.worst.4AD)

SA_pred1 <- dplyr::left_join(SA_data1, pred_data1, by = c("Title" = "Tetra"))
SA_pred2 <- dplyr::left_join(SA_data2, pred_data2, by = c("Title" = "Tetra"))
SA_pred3 <- dplyr::left_join(SA_data3, pred_data3, by = c("Title" = "Tetra"))
SA_pred4 <- dplyr::left_join(SA_data4, pred_data4, by = c("Title" = "Tetra"))

SA_pred1AD <- dplyr::left_join(SA_data1AD, pred_data1AD, by = c("Title" = "Tetra"))
SA_pred2AD <- dplyr::left_join(SA_data2AD, pred_data2AD, by = c("Title" = "Tetra"))
SA_pred3AD <- dplyr::left_join(SA_data3AD, pred_data3AD, by = c("Title" = "Tetra"))
SA_pred4AD <- dplyr::left_join(SA_data4AD, pred_data4AD, by = c("Title" = "Tetra"))

SA_pred1 <- SA_pred1[SA_pred1$Title %in% rownames(data_seq2)[1:(nrow(pred_data1)*0.05)],]
SA_pred2 <- SA_pred2[SA_pred2$Title %in% rownames(data_seq2)[1:(nrow(pred_data2)*0.05)],]
SA_pred3 <- SA_pred3[SA_pred3$Title %in% rownames(data_seq3)[1:(nrow(pred_data3)*0.05)],]
SA_pred4 <- SA_pred4[SA_pred4$Title %in% rownames(data_seq4)[1:(nrow(pred_data4)*0.05)],]

SA_pred1AD <- SA_pred1AD[SA_pred1AD$Title %in% rownames(data_seq1AD)[1:(nrow(pred_data1AD)*0.05)],]
SA_pred2AD <- SA_pred2AD[SA_pred2AD$Title %in% rownames(data_seq2AD)[1:(nrow(pred_data2AD)*0.05)],]
SA_pred3AD <- SA_pred3AD[SA_pred3AD$Title %in% rownames(data_seq3AD)[1:(nrow(pred_data3AD)*0.05)],]
SA_pred4AD <- SA_pred4AD[SA_pred4AD$Title %in% rownames(data_seq4AD)[1:(nrow(pred_data4AD)*0.05)],]

# PLOTS

library(ggplot2)
library(dplyr)

calc_metrics <- function(selection, data, best_data) {
  
  test_dock <- data[order(data$Pred, decreasing = T),]
  test_dock <- test_dock[1:selection,]
  test_best <- best_data
  
  test_coocur <- test_dock$Tetra %in% test_best$Title
  test_coocur <- test_dock[test_coocur,]
  
  p1 <- nrow(test_coocur)/nrow(test_best)
  p2 <- nrow(test_coocur)/nrow(test_dock)
  
  return(list(p1=p1, p2=p2))
}

calc_metrics_list <- function(datasets, best_datasets, selection) {
  
  results_list <- mapply(function(data, best_data) {
    
    results <- lapply(selection, function(sele) calc_metrics(sele, data, best_data))
    results_df <- data.frame(matrix(unlist(results), nrow=length(results), byrow=T))
    names(results_df) <- c("p1", "p2")
    results_df$selection <- selection
    return(results_df)
    
  }, datasets, best_datasets, SIMPLIFY = FALSE)
}


datasets1 <- list(   pred_data1, pred_data2, pred_data3, pred_data4)
best_datasets1 <- list(SA_pred1,   SA_pred2,   SA_pred3,   SA_pred4)

datasets2 <- list(   pred_data1AD, pred_data2AD, pred_data3AD, pred_data4AD)
best_datasets2 <- list(SA_pred1AD,   SA_pred2AD,   SA_pred3AD,   SA_pred4AD)

selection <- seq(0,50000, 500)

results_df1 <- calc_metrics_list(datasets1, best_datasets1, selection)
names(results_df1) <- c("CHIKV", "DENV", "WNV", "ZIKV")

results_df2 <- calc_metrics_list(datasets2, best_datasets2, selection)
names(results_df2) <- c("CHIKV_AD", "DENV_AD", "WNV_AD", "ZIKV_AD")

average_results1 <- Reduce("+", results_df1) / length(results_df1)
average_results1$group <- "Openeye"

average_results2 <- Reduce("+", results_df2) / length(results_df2)
average_results2$group <- "AutoDock"

final_results_df <- rbind(average_results1, average_results2)

q1 <- ggplot(final_results_df, aes(x=selection, y=p1, color=group)) + 
  geom_line(linewidth = 3) +
  xlab("ML Selected Peptides") +
  ylab("Concurrence")+
  scale_y_continuous(labels = scales::percent)+
  theme_bw()+
  theme(text = element_text(size = 50, color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=4.5),
        axis.text = element_text(size = 50, color = "black"),
        panel.grid = element_line(color = "gray"),
        legend.title=element_blank())

q1

ggsave(plot = q1, device = "png", filename ="./Not_in_paper/Concurrence.png",width = 20, height = 10, dpi = 300)

write.csv(final_results_df, "./Not_in_paper/Structural_Analysis.csv", row.names = F)

# Table 5
table5 <- final_results_df[final_results_df$selection == 100|
                             final_results_df$selection == 500|
                             final_results_df$selection == 1000|
                             final_results_df$selection == 2000|
                             final_results_df$selection == 4000|
                             final_results_df$selection == 5000|
                             final_results_df$selection == 6000|
                             final_results_df$selection == 7000|
                             final_results_df$selection == 8000|
                             final_results_df$selection == 16000|
                             final_results_df$selection == 32000|
                             final_results_df$selection == 50000,]


write.csv(table5, "./Results/Table_5.csv", row.names = F)

# Plot Figure S5

data_preds <- list(pred_data1, pred_data2, pred_data3, pred_data4, 
                   pred_data1AD, pred_data2AD, pred_data3AD, pred_data4AD)
datasets <- c("CHIKV",    "DENV",    "WNV",    "ZIKV",
              "CHIKV_AD", "DENV_AD", "WNV_AD", "ZIKV_AD")

peptide_numbers <- c(100, 1000, 2000, 4000, 6000, 8000, 16000, 32000, 50000)

for (i in 1:length(data_preds)) {
  
  res <- data_preds[[i]]
  ds <- datasets[i]
  

  data_subsets <- list()
  data_pred <- res[order(res$Pred, decreasing = T),]
  
  # Loop over the peptide numbers
  for (j in 1:length(peptide_numbers)) {
    # Subset the data frame and store it in the list
    data_subsets[[j]] <- data_pred[1:peptide_numbers[j], ]
    
    # Create a new column in this subset to store the peptide number (for facetting the plot)
    data_subsets[[j]]$peptide_number <- peptide_numbers[j]
  }
  
  # Combine all the subsets into one data frame
  all_data <- do.call(rbind, data_subsets)
  
  # Plot the boxplots
  p <- ggplot(all_data, aes(x=factor(peptide_number), y=score)) +
    geom_boxplot() +
    stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 2)), vjust = -0.5, size = 7)+
    annotate("text", x = "6000", y = max(all_data$score), label = ds, size = 40, hjust = 0.5, vjust = 1)+
    labs(x="N. Selected peptides", y="Score") +
    scale_x_discrete(limits = rev(levels(factor(all_data$peptide_number)))) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          text = element_text(size = 36),
          axis.title = element_text(size = 36), 
          axis.text = element_text(size = 36, color = "black"),
          panel.border = element_rect(linewidth = 3),
          panel.grid = element_line(color = "gray"))  
  
  
  v <- letters[i]
  print(p)
  ggsave(plot = p, device = "png", paste0("./Supplementary/Figure_S5.", v, ".png"),width = 20, height = 17.5, dpi = 600)
}

