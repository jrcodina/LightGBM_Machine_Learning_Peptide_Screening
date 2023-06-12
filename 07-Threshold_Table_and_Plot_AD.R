library(lightgbm)
library(dplyr)
library(caret)

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

model1 <- readRDS.lgb.Booster("./Results/LightGBM_MODELS/model_AD_CHIKV.RDs")
model2 <- readRDS.lgb.Booster("./Results/LightGBM_MODELS/model_AD_DENV.RDs")
model3 <- readRDS.lgb.Booster("./Results/LightGBM_MODELS/model_AD_WNV.RDs")
model4 <- readRDS.lgb.Booster("./Results/LightGBM_MODELS/model_AD_ZIKV.RDs")

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
data_seq3 <- prepare.datasets_HO(data_sequence2)
data_seq4 <- prepare.datasets_HO(data_sequence3)
data_seq5 <- prepare.datasets_HO(data_sequence4)

extract_cmTable <- function(data_seq, reference_data, model){
  
  train_idx <- sample(1:nrow(data_seq), 0.01 * nrow(data_seq))
  
  train_data <- data_seq[train_idx, ]
  test_data <- data_seq[-train_idx, ]
  
  X_train <- as.matrix(train_data[, !(names(train_data) %in% c("TARGET", "Tetra"))]) 
  y_train <- train_data$TARGET
  
  X_test <- test_data[, !(names(test_data) %in% c("TARGET", "Tetra"))]
  y_test <- test_data$TARGET
  
  pred <- predict(model, as.matrix(X_test), reshape = T)
  
  real_vs_pred <- data.frame(Tetra = rownames(test_data),  Real = y_test, Pred = pred)
  real_vs_pred <- dplyr::left_join(real_vs_pred, reference_data[,1:2])
  
  thresholds <- seq(0.05, .95, 0.05)
  predictions_th <- list()
  
  for (i in 1:length(thresholds)) {
    tmp <- ifelse(real_vs_pred$Pred >= thresholds[i], 1, 0)
    predictions_th[[i]] <- tmp
    names(predictions_th)[i] <- paste0("Threshold_", thresholds[i])
  }
  predictions_th <- data.frame(do.call(cbind, predictions_th))
  
  real_vs_pred <- cbind(real_vs_pred, predictions_th)
  
  cmLGBMs <- list()
  
  for (i in 1:length(thresholds)) {
    cmLGBM <- caret::confusionMatrix(reference = as.factor(real_vs_pred$Real), as.factor(real_vs_pred[,4+i]), positive = "1")
    cmLGBMs[[i]] <- cmLGBM
  }
  names(cmLGBMs) <- thresholds
  
  cmTable <- c()
  for (i in 1:length(cmLGBMs)) {
    cmTable[[i]] <- data.frame(Threshold = thresholds[i],
                               Accuracy = cmLGBMs[[i]][["overall"]][["Accuracy"]],
                               PPV = cmLGBMs[[i]][["byClass"]][["Pos Pred Value"]],
                               NPV = cmLGBMs[[i]][["byClass"]][["Neg Pred Value"]],
                               TP = cmLGBMs[[i]][["table"]]["1","1"],
                               TN = cmLGBMs[[i]][["table"]]["0","0"],
                               FP = cmLGBMs[[i]][["table"]]["1","0"],
                               FN = cmLGBMs[[i]][["table"]]["0","1"],
                               Sensitivity = cmLGBMs[[i]][["byClass"]][["Sensitivity"]],
                               Specificity = cmLGBMs[[i]][["byClass"]][["Specificity"]],
                               Pred.as.Positive = sum(cmLGBMs[[i]][["table"]]["1",]),
                               Pred.as.Negative = sum(cmLGBMs[[i]][["table"]]["0",]),
                               F1 = cmLGBMs[[i]][["byClass"]][["F1"]]
    )
  }
  
  cmTable <- as.data.frame(do.call(rbind, cmTable))
  cmTable[is.na(cmTable)] <- 0
  return(cmTable)
}

cmCHIKV <- extract_cmTable(data_seq = data_seq1, reference_data = data1, model = model1)
cmDENV  <- extract_cmTable(data_seq = data_seq2, reference_data = data2, model = model2)
cmWNV   <- extract_cmTable(data_seq = data_seq3, reference_data = data3, model = model3)
cmZIKV  <- extract_cmTable(data_seq = data_seq4, reference_data = data4, model = model4)

write.csv(cmZIKV, "./Results/ConfussionMatrix/cmZIKV_AD.csv", row.names = F)
write.csv(cmDENV, "./Results/ConfussionMatrix/cmDENV_AD.csv", row.names = F)
write.csv(cmCHIKV, "./Results/ConfussionMatrix/cmCHIKV_AD.csv", row.names = F)
write.csv(cmWNV, "./Results/ConfussionMatrix/cmWNV_AD.csv", row.names = F)

# Plots

cmTables <- list(cmCHIKV, cmDENV, cmWNV, cmZIKV)
datasets <- c("AD_CHIKV", "AD_DENV", "AD_WNV", "AD_ZIKV")

for (i in 1:length(cmTables)) {
  dataset <- datasets[i]
  cm <- cmTables[[i]]
  
  p <- ggplot(data = cm, aes(x = Threshold))+
    geom_col(aes(y = TP + FP, fill = "True Positives"), color = "black") +
    geom_col(aes(y = FP, fill = "False Positives"), color = "black") +
    geom_line(aes(y = PPV*max(Pred.as.Positive), color = "PPV"), linewidth = 3)+
    geom_point(aes(y = PPV*max(Pred.as.Positive), color = "PPV"), size = 5)+
    geom_text(aes(y = TP + FP, label = TP+FP, fontface = "bold"), vjust = -0.3, hjust = 0.35, size = 10)+
    geom_text(aes(x = 0.5, y = max(TP + FP), label = dataset), hjust = 0.5, vjust =1, size = 40)+
    scale_y_continuous(breaks = seq(0, 125000, 25000), labels = scales::comma, name = "Predicted as Positive", 
                       sec.axis = sec_axis(~./max(cm$Pred.as.Positive), name = "PPV", breaks = seq(0, 1, 0.1)))+
    scale_fill_manual(values = c("True Positives" = "steelblue", "False Positives" = "orange"),
                      guide = guide_legend(title = ""))+
    scale_color_manual(values = c("PPV" = "red"), 
                       name = "", 
                       labels = c("PPV"))+
    scale_x_continuous(breaks = seq(0, 1, 0.1))+
    theme_bw()+
    theme(text = element_text(size = 45),
          axis.title = element_text(size = 45), 
          axis.text = element_text(size = 45, color = "black"),
          panel.border = element_rect(linewidth = 4),
          panel.grid = element_line(color = "gray"),
          axis.line.y.right = element_line(color = "red", linewidth = 3), 
          axis.ticks.y.right = element_line(color = "red"),
          axis.text.y.right = element_text(color = "red"),
          axis.title.y.right = element_text(color = "red"))+
    guides(fill = "none", color = "none")
  
  v <- letters[i]
  ggsave(plot = p, device = "png", paste0("./Results/Figures/Figure_4.2.",v, ".png"),width = 20, height = 17.5, dpi = 600)
}

p <- ggplot(data = cm, aes(x = Threshold))+
  geom_col(aes(y = TP + FP, fill = "True Positives"), color = "black") +
  geom_col(aes(y = FP, fill = "False Positives"), color = "black") +
  geom_line(aes(y = PPV*max(Pred.as.Positive), color = "PPV"), linewidth = 3)+
  geom_point(aes(y = PPV*max(Pred.as.Positive), color = "PPV"), size = 5)+
  geom_text(aes(y = TP + FP, label = TP + FP), vjust = -0.5, size = 6) +
  scale_y_continuous(name = "Predicted as Positive", 
                     sec.axis = sec_axis(~./max(cm$Pred.as.Positive), name = "PPV", breaks = seq(0, 1, 0.1)))+
  scale_fill_manual(values = c("True Positives" = "steelblue", "False Positives" = "orange"),
                    guide = guide_legend(title = ""))+
  scale_color_manual(values = c("PPV" = "red"), 
                     name = "", 
                     labels = c("PPV"))+
  scale_x_continuous(breaks = seq(0, 1, 0.1))+
  labs(title = dataset)+
  theme_bw()+
  theme(text = element_text(size = 36),
        axis.title = element_text(size = 36), 
        axis.text = element_text(size = 36, color = "black"),
        panel.border = element_rect(linewidth = 3),
        panel.grid = element_line(color = "gray"),
        axis.line.y.right = element_line(color = "red"), 
        axis.ticks.y.right = element_line(color = "red"),
        axis.text.y.right = element_text(color = "red"),
        axis.title.y.right = element_text(color = "red"))
legend <- cowplot::get_legend(p)
ggsave(plot = legend, device = "png", "./Results/Figures/Figure_4_legend.png",width = 10, height = 10, dpi = 1200)
