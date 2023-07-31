#dt = 17.25 mins
t0 <- Sys.time()

library(lightgbm)
library(dplyr)
library(tidyr)
library(data.table)

results <- list(ZIKV = list(readRDS.lgb.Booster("./Results/LightGBM_MODELS/model_ZIKV.RDs")),
                DENV = list(readRDS.lgb.Booster("./Results/LightGBM_MODELS/model_DENV.RDs")),
                CHIKV = list(readRDS.lgb.Booster("./Results/LightGBM_MODELS/model_CHIKV.RDs")),
                WNV = list(readRDS.lgb.Booster("./Results/LightGBM_MODELS/model_WNV.RDs")))

# Initialize a list to store the average feature importance for each dataset
avg_feature_importance <- vector("list", 4)

for (i in 1:4) {
  # Extract the feature importance of each model in the current dataset
  feature_importance_list <- lapply(results[[i]], function(model) {
    return(lightgbm::lgb.importance(model))
  })
  
  # Combine the feature importance data.frames into a single data.frame
  combined_importance <- do.call(rbind, feature_importance_list)
  
  # Calculate the average feature importance per feature
  avg_importance <- aggregate(combined_importance[, c("Gain", "Cover", "Frequency")],
                              by = list(Feature = combined_importance$Feature),
                              FUN = mean)
  
  # Store the average feature importance data.frame in the list
  avg_feature_importance[[i]] <- avg_importance
  rm(avg_importance, combined_importance)
}

names(avg_feature_importance) <- c("CHIKV", "DENV", "WNV", "ZIKV")
# 'avg_feature_importance' now contains the average feature importance for each dataset

# Plots
library(ggplot2)

create_ggplot <- function(avg_importance, dataset_name) {
  plot <- ggplot(avg_importance, aes(x = reorder(Feature, -Gain), y = Gain)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = paste("Feature Importance for ", dataset_name, " Dataset"),
         x = "Feature",
         y = "Average Gain")
  
  return(plot)
}

# Plot the average feature importance for each dataset
for (i in 1:4) {
  dataset_name <- names(avg_feature_importance)[i]
  plot <- create_ggplot(avg_feature_importance[[i]], dataset_name)
  print(plot)
}


# Combine the average feature importances into a single data frame
all_avg_importance <- do.call(rbind, lapply(1:4, function(i) {
  df <- avg_feature_importance[[i]]
  df$Dataset <- names(avg_feature_importance)[i]
  return(df)
}))

# Create a ggplot object for the combined average feature importances
plot_all <- ggplot(all_avg_importance, aes(x = reorder(Feature, -Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Average Feature Importance for All Datasets",
       x = "Feature",
       y = "Average Gain") +
  facet_wrap(~ Dataset, scales = "free_y", ncol = 2)

plot_all

ggsave("./Supplementary/Figure_S4.b.png", plot_all, width = 25, height = 6.66, dpi = 600)
write.csv(all_avg_importance, "./Supplementary/Table_S3.2.csv")


# Top 5 important features
top_features_df <- data.frame()

# Loop through each dataset and extract top 5 features
for(i in 1:length(avg_feature_importance)) {
  dataset <- avg_feature_importance[[i]]
  top_features <- dataset %>%
    arrange(desc(Gain)) %>%  
    head(5)                   # Get top 5
  
  top_features$Dataset <- names(avg_feature_importance)[i]
  
  top_features_df <- rbind(top_features_df, top_features)
}

top_features_df <- top_features_df[, c("Dataset", "Feature", "Gain")]

combined_df <- do.call(rbind, avg_feature_importance)

combined_df$Dataset <- rep(names(avg_feature_importance), sapply(avg_feature_importance, nrow))

feature_stats <- combined_df %>%
  group_by(Feature) %>%
  summarise(Avg_Gain = mean(Gain), SD_Gain = sd(Gain), .groups = "drop")

# Get the top 5 features based on average Gain
top_features_all <- feature_stats %>%
  arrange(desc(Avg_Gain)) %>%
  head(5)

write.csv(top_features_df, "./Supplementary/Table_4.2.csv", row.names = F)
write.csv(top_features_all, "./Not_in_paper/Top_Features_Average_AD.csv", row.names = F)

