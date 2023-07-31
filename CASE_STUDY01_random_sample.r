library(dplyr)
library(lightgbm)

tetra <- expand.grid(rep(list(Peptides::aaList()),4))

tetra <- do.call(paste0, tetra[, 1:4])

set.seed(1234)
tetra <- tetra[sample(1:160000, 1600)]

tetra_sample <- paste0(tolower(tetra_sample), "g") # Because lower case means extanded conformation, while uppercase means a-helical conformation. The alpha fold test that we have done represents the pentapeptides as extended, so we are doing extended.

cat(tetra_sample, file = "tetra_G_sample.txt", sep = "\n", append = FALSE)
