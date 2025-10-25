rm(list=ls())
# library()
parse_pt_file <- function(filename, print_res=TRUE) {
  # Read all lines from the file
  lines <- readLines(filename)
  
  # Initialize empty vectors
  all_temps <- numeric()
  all_iterations <- numeric()
  all_swap_rate <- numeric()
  all_seconds <- numeric()
  
  # Process each block
  i <- 1
  while (i <= length(lines)) {
    # Skip empty lines
    if (lines[i] == "") {
      i <- i + 1
      next
    }
    
    # Check if this is the start of a new block
    if (grepl("^alg:", lines[i])) {
      # Extract data from this block
      block_data <- list()
      
      # Process each line in the block
      while (i <= length(lines) && lines[i] != "") {
        if(grepl("^temp_ini:", lines[i])){
          
        }
        else if (grepl("^\\s*temps:", lines[i])) {
          temp_str <- gsub("^\\s*temps:\\s*", "", lines[i])
          temp_str <- gsub(",$", "", temp_str)  # Remove trailing comma
          block_data$temps <- as.numeric(strsplit(temp_str, ",")[[1]])
        }
        else if (grepl("^\\s*iterations:", lines[i])) {
          iter_str <- gsub("^\\s*iterations:\\s*", "", lines[i])
          iter_str <- gsub(",$", "", iter_str)  # Remove trailing comma
          block_data$iterations <- as.numeric(strsplit(iter_str, ",")[[1]])
        }
        else if (grepl("^\\s*swap_rate:", lines[i])) {
          swap_str <- gsub("^\\s*swap_rate:\\s*", "", lines[i])
          swap_str <- gsub(",$", "", swap_str)  # Remove trailing comma
          block_data$swap_rate <- as.numeric(strsplit(swap_str, ",")[[1]])
        }
        else if (grepl("^seconds:", lines[i])) {
          block_data$seconds <- as.numeric(gsub("^seconds:\\s*", "", lines[i]))
        }
        i <- i + 1
      }
      
      # Add data to main vectors, handling duplicates for temps and iterations
      if (!is.null(block_data$temps)) {
        if (length(all_temps) == 0) {
          # First block - include all temps
          all_temps <- c(all_temps, block_data$temps)
        } else {
          # Subsequent blocks - skip first temp (duplicate of last from previous block)
          all_temps <- c(all_temps, block_data$temps[-1])
        }
      }
      
      if (!is.null(block_data$iterations)) {
        if (length(all_iterations) == 0) {
          # First block - include all iterations
          all_iterations <- c(all_iterations, block_data$iterations)
        } else {
          # Subsequent blocks - skip first iteration (duplicate of last from previous block)
          all_iterations <- c(all_iterations, block_data$iterations[-1])
        }
      }
      
      # For swap_rate and seconds, include all data
      if (!is.null(block_data$swap_rate)) {
        all_swap_rate <- c(all_swap_rate, block_data$swap_rate)
      }
      
      if (!is.null(block_data$seconds)) {
        all_seconds <- c(all_seconds, block_data$seconds)
      }
    } else {
      i <- i + 1
    }
  }
  
  if(print_res){
    cat("Temps:\n");
    print(paste0(all_temps,collapse = ","));
    cat("Iter: \n");
    print(paste0(all_iterations,collapse = ","));
    cat("Swap rate:\n");
    print(paste0(all_swap_rate,collapse = ","));
  }
  # Return as a list of vectors
  return(list(
    temps = all_temps,
    iterations = all_iterations,
    swap_rate = all_swap_rate,
    seconds = all_seconds
  ))
}

id <- 182
dat <- parse_pt_file(file.path(getwd(),"results",paste0("temperatures_id_",id,".txt")))
            