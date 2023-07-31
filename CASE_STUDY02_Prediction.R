# Load necessary package
library(stringr)

# Define filenames
filenames <- paste0("WNV_", 1:10, ".out")

# Initialize a list to hold best affinities from all files
all_best_affinities <- list()

# Loop over all files
for (filename in filenames) {

  # Read in file
  file_text <- readLines(filename)

  # Find lines where affinity data starts
  affinity_starts <- grep("   1       ", file_text)

  # Calculate number of peptides
  num_peptides <- length(affinity_starts)

  # Initialize vector to hold best affinities
  best_affinities <- numeric(num_peptides)

  # Loop over all peptides
  for (i in 1:num_peptides) {
    # Extract line with best affinity data for this peptide
    peptide_data <- file_text[affinity_starts[i]]

    # Extract affinity (assuming it always starts at 11th character and end at 17th character)
    affinity <- as.numeric(str_sub(peptide_data, 11, 17))

    # Save the affinity
    best_affinities[i] <- affinity
  }

  # Create a data frame from the best affinities
  best_affinities_df <- data.frame(
    "File" = filename,
    "Peptide" = 1:num_peptides,
    "BestAffinity" = best_affinities
  )

  # Append to the overall data frame
  all_best_affinities <- bind_rows(all_best_affinities, best_affinities_df)
}

# Print best affinities from all files
print(all_best_affinities)

tetra_sample <- read.csv("./tetra_G_sample.txt", header = F)

results_sample <- cbind(tetra_sample, all_best_affinities)

# LightGBM prediction

library(lightgbm)
Best_Params <- read.csv("./Results/Best_params.csv")
results_sample <- cbind(tetra_sample, all_best_affinities)

data <- results_sample[,c(1,4)]
names(data) <- c("Sequence", "X")

data <- data[order(data$X),]

data_prop <- PCAprop(custom.list = T, PeList = data$Sequence, norm = T)
data_prop <- dplyr::left_join(data, data_prop)

prepare.datasets <- function(x, positive_negative = c(0.2, 0.8)){
  x <- x[x$Sequence != "",]
  x <- na.omit(x)
  x <- x[order(x[,2]),]

  xn <- x[,-c(1:2)]
  rownames(xn) <- x$Sequence

  x1 <- replicate(nrow(xn)*positive_negative[1],1)
  x2 <- replicate(nrow(xn)*positive_negative[2],0)

  y <- data.frame(y = c(x1, x2))
  y <- rbind(y, data.frame(y = rep(0, nrow(xn) - nrow(y))))

  XYn <- as.data.frame(cbind(y, xn))

  names(XYn)[1] <- "TARGET"
  return(XYn)
}

data_seq <- prepare.datasets(data_prop)


train_idx <- sample(1:nrow(data_seq), 0.75 * nrow(data_seq))

train_data <- data_seq[train_idx, ]
test_data <- data_seq[-train_idx, ]

X_train <- as.matrix(train_data[, !(names(train_data) %in% c("TARGET", "Sequence"))])
y_train <- train_data$TARGET

X_test <- test_data[, !(names(test_data) %in% c("TARGET", "Sequence"))]
y_test <- test_data$TARGET

# Set up the parameters
params <- list(
  objective = "binary",
  metric = "auc",
  boosting = "gbdt",
  scale_pos_weight = Best_Params$scale_pos_weight,
  learning_rate = Best_Params$learning_rate,
  num_leaves = Best_Params$num_leaves,
  feature_fraction = Best_Params$feature_fraction,
  pos_bagging_fraction = Best_Params$pos_bagging_fraction,
  neg_bagging_fraction = Best_Params$neg_bagging_fraction,
  bagging_freq = Best_Params$bagging_freq,
  min_data_in_leaf = Best_Params$min_data_in_leaf,
  max_depth = Best_Params$max_depth

)

dtrain = lgb.Dataset(data = as.matrix(X_train), label = y_train)
valids = list(validation = lgb.Dataset(data = as.matrix(X_test), label = y_test))

# Train the model
model <- lgb.train(
  data = dtrain,
  valids = valids,
  params = params,
  nrounds = 5000,
  early_stopping_rounds = 50)

to_pred <- PCAprop(custom.list = T, PeList = OE_res$X3N40OnlyErec_1, norm = T)
X_topred <- to_pred$Sequence
Y_topred <- as.matrix(to_pred[,-1])

pred <- predict(model, as.matrix(Y_topred), reshape = T)

prediction <- cbind(X_topred, pred)
prediction <- prediction[order(prediction[,2], decreasing = T),]

selection_5 <- prediction[1:8000,]
selection_worst_5 <- prediction[160000:152000,]

write.csv(prediction, "./Table_S5.csv", row.names = F)
write.csv(results_sample, "./Table_S4.csv", row.names = F)
lightgbm::saveRDS.lgb.Booster(model, "./Model_Tetra_G.lgb")
write.csv(selection_5, "./20230706_Selected_8000.csv")
cat(selection_5[,1], file = "tetra_G_selected.txt", sep = "\n", append = FALSE)
cat(selection_worst_5[,1], file = "tetra_G_selected_worst.txt", sep = "\n", append = FALSE)


# assuming 'pred' is your vector of values
h <- hist(pred, plot=FALSE, breaks = seq(0,1,.1))  # create histogram data

# plot histogram
plot(h, xlab="pred", main="Histogram with frequencies")

# add frequencies on top of the bars
text(x=h$mids, y=h$counts, labels=h$counts, pos=3, cex=0.8, col="red")
