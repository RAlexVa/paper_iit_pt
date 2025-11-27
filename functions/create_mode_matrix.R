create_binary_matrix <- function(n, probabilities=c(0.1, 0.3, 0.6, 0.8)) {
  # Create empty matrix
  mat <- matrix(0, nrow = n, ncol = n)
  
  # Define number of groups
  groups <- length(probabilities)
  # Calculate column indices for each chunk
  chunk <- n %/% groups # Divide columns in groups of roughly the same size
  
  col_indices <- list()
  for(i in 1:(groups-1)){
    col_indices[[i]] <- (1 + chunk*(i-1)):(chunk*i)
  }
  
  col_indices[[groups]] <- (1 + chunk*(groups-1)):(n)

  
  
  # Fill each quarter with random 0s and 1s
  for (i in 1:groups) {
    cols <- col_indices[[i]]
    p <- probabilities[i]
    # Generate random matrix for this quarter using binomial distribution
    mat[, cols] <- matrix(rbinom(n * length(cols), 1, p), nrow = n)
  }
  
  return(mat)
}

write_matrix_to_file_advanced <- function(mat, filename, format = "space") {
  # Ensure the filename has .txt extension
  if (!grepl("\\.txt$", filename)) {
    filename <- paste0(filename, ".txt")
  }
  
  # Choose separator based on format
  separator <- switch(format,
                      "space" = " ",
                      "comma" = ",",
                      "tab" = "\t",
                      " ")
  
  # Write matrix to file
  write.table(mat, file = filename, 
              row.names = FALSE, col.names = FALSE,
              sep = separator)
  
  cat("Matrix written to:", filename, "\n")
  cat("Matrix dimensions:", nrow(mat), "x", ncol(mat), "\n")
  cat("File format:", format, "separated\n")
}


compute_column_differences <- function(mat) {
  n_cols <- ncol(mat)
  
  # Initialize results
  results <- list()
  combination_count <- 0
  
  # Calculate for all unique pairs of columns
  for (i in 1:(n_cols - 1)) {
    for (j in (i + 1):n_cols) {
      diff_sum <- sum(abs(mat[, i] - mat[, j]))
      combination_count <- combination_count + 1
      
      results[[combination_count]] <- list(
        col1 = i,
        col2 = j,
        diff_sum = diff_sum
      )
    }
  }
  
  # Convert to data frame for easier viewing
  results_df <- do.call(rbind, lapply(results, as.data.frame))
  
  return(results_df)
}




# set.seed(1)
# a <- create_binary_matrix(100)
# set.seed(1)
# b <- create_binary_matrix(100)
identical(a,b)


for(i in 1:100){
  for(dim in c(500,800,1000)){
    print(paste(c("seed:",i,"dim:",dim),collapse = " "))
    set.seed(i);
    write_matrix_to_file_advanced(
      create_binary_matrix(dim),
      filename=paste0("inputs/mode_matrix/","mat_",dim,"_",i,".txt")
    )
  }
}

Rcpp::sourceCpp("functions/highdim_parallel_sep_multimodal.cpp")

test <- readMatrix("inputs/mode_matrix/mat_500_3.txt")

dist <- compute_column_differences(test)
