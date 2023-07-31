# Load necessary package
library(dplyr)

# Set the path to the directory containing the files
path <- "./results_ADCP/bad"

# Get all file names
file_list <- list.files(path, pattern = ".out")

# Create an empty dataframe to store results
results_df <- data.frame(filename = character(), molecule_number = numeric(), affinity = numeric(), error_type = character(), time = numeric())

# Loop over all files
for(file in file_list){
  # Read the file
  text <- readLines(paste0(path, "/", file))

  # Add a line at the end to ensure the last molecule is processed correctly
  text <- c(text, "Detected 16 cores")

  # Split the file into different molecules
  split_idx <- grep("Detected 16 cores", text)

  while (length(split_idx) < 161) {
    text <- c(text, "Detected 16 cores")
    split_idx <- grep("Detected 16 cores", text)
  }

  molecules <- split(text, cut(seq_along(text), split_idx, labels = FALSE))

  # Loop over all molecules in the file
  for(i in 1:length(molecules)){
    molecule <- molecules[[i]]

    # Count the number of lines
    n_lines <- length(molecule)

    # Check for errors
    if (n_lines == 1) {
      error <- "incomplete output"
    } else {
      if(!any(grepl("   1       ", molecule))){
        # Check for "no output" error
        if(!grepl("clean up unzipped map files", molecule[n_lines -1]) |
           !grepl("", molecule[n_lines -1])){
          error <- "incomplete output"
        }
        # Check for "file already exists" error
        else if(any(grepl("ERROR: output file exists", molecule))){
          error <- "file already exists"
        }
        else{
          error <- NA
        }
      }
      else{
        error <- NA
      }
    }

    if (any(grepl("   1       ", molecule))) {
      affinity <- molecule[grepl("   1       ", molecule)]
      affinity <- as.numeric(substr(affinity, 11, 17))
    }else{
      affinity <- NA
    }

    if (any(grepl("Docking performed in ", molecule))){
      time <- molecule[grepl("Docking performed in ", molecule)]
      time <- as.numeric(substr(time, 21, 28))
    }else{
      time <- NA}

    # Add the information to the dataframe
    results_df <- results_df %>%
      add_row(filename = file, molecule_number = i, affinity = affinity, error_type = error, time = time)
  }
}

# Print the results
print(results_df)

print(paste("Errors:", length(which(!is.na(results_df$error_type)))))

all_best_affinities <- results_df

write.csv(all_best_affinities, row.names = F, "./20230717_Results_ADCP_Worse.csv")

sequences <- read.csv("./tetra_G_selected_worst.txt", header = F)
all_best_affinities <- cbind(Sequence = sequences[1:8000,], all_best_affinities)

write.csv(all_best_affinities, row.names = F, "./Table_S7.csv")
