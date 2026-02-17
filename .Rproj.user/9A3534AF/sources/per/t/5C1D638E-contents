rm(list=ls())
# library()

########## Parse a single file #####
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
id <- 265
dat <- parse_pt_file(file.path(getwd(),"results",paste0("temperatures_id_",id,".txt")))
            

##### Parse multiple files #####
### To concatenate many files in a single output matrix
parse_pt_file <- function(filename, print_res=FALSE) {
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
  
  # Return as a data frame with equal length vectors
  max_length <- max(length(all_temps), length(all_iterations), 
                    length(all_swap_rate), length(all_seconds))
  
  # Pad shorter vectors with NA to make them all the same length
  temps_padded <- c(all_temps, rep(NA, max_length - length(all_temps)))
  iterations_padded <- c(all_iterations, rep(NA, max_length - length(all_iterations)))
  swap_rate_padded <- c(all_swap_rate, rep(NA, max_length - length(all_swap_rate)))
  seconds_padded <- c(all_seconds, rep(NA, max_length - length(all_seconds)))
  
  return(data.frame(
    temps = temps_padded,
    iterations = iterations_padded,
    swap_rate = swap_rate_padded,
    seconds = seconds_padded
  ))
}

# Function to process multiple files and combine them
process_multiple_files <- function(ids, vector_names = NULL) {
  # If vector_names not provided, use default names
  if (is.null(vector_names)) {
    vector_names <- paste0("vector_", seq_along(ids))
  }
  
  # Check that ids and vector_names have the same length
  if (length(ids) != length(vector_names)) {
    stop("ids and vector_names must have the same length")
  }
  
  # Initialize list to store results
  results_list <- list()
  
  # Process each file
  for (j in seq_along(ids)) {
    id <- ids[j]
    vector_name <- vector_names[j]
    
    filename <- file.path(getwd(), "results", paste0("temperatures_id_", id, ".txt"))
    
    # Check if file exists
    if (!file.exists(filename)) {
      warning(paste("File not found:", filename))
      next
    }
    
    # Parse the file
    dat <- parse_pt_file(filename, print_res = FALSE)
    
    # Add identification columns
    dat$id <- id
    dat$vector_name <- vector_name
    
    # Add to results list
    results_list[[j]] <- dat
  }
  
  # Combine all results using rbind
  final_matrix <- do.call(rbind, results_list)
  
  return(final_matrix)
}

temps_table <- process_multiple_files(233:264)

write.table(temps_table, 
            file = "results/output_results.txt",
            sep = "\t",           # tab-separated
            row.names = FALSE,    # don't include row names
            col.names = TRUE,     # include column names
            quote = FALSE,        # don't put quotes around strings
            na = "NA")            # how to represent missing values



##### Process multiple files and transpose the results ##### 
temps_table <- process_multiple_files(299:302)
# First, let's create a function to process each ID's data
process_single_id <- function(result_matrix, target_id) {
  # Filter rows for this specific ID and drop seconds and vector_name columns
  id_data <- result_matrix[result_matrix$id == target_id, c("id", "temps", "iterations", "swap_rate")]
  
  # Remove rows with all NAs (from padding)
  id_data <- id_data[!is.na(id_data$temps), ]
  
  # Create the transposed matrix
  # Rows will be: temps, iterations, swap_rate
  # Columns will be the temperature values
  
  # Get unique temperature values (they should be the same across all metrics for this ID)
  temp_values <- unique(id_data$temps)
  temp_values <- temp_values[!is.na(temp_values)]
  
  # Create the transposed matrix
  transposed_data <- data.frame(
    id = target_id,
    metric = c("temps", "iterations", "swap_rate")
  )
  
  # Add columns for each temperature value
  for (i in seq_along(temp_values)) {
    temp_col <- numeric(3)
    
    # Get values for this temperature
    temp_val <- temp_values[i]
    
    # temps row (should be the temperature value itself)
    temp_col[1] <- temp_val
    
    # iterations row
    iter_val <- id_data$iterations[id_data$temps == temp_val]
    temp_col[2] <- ifelse(length(iter_val) > 0, iter_val[1], NA)
    
    # swap_rate row
    swap_val <- id_data$swap_rate[id_data$temps == temp_val]
    temp_col[3] <- ifelse(length(swap_val) > 0, swap_val[1], NA)
    
    # Add column to transposed data
    col_name <- paste0("temp_", i)
    transposed_data[[col_name]] <- temp_col
  }
  
  return(transposed_data)
}

# Function to process all IDs and combine them
process_all_ids <- function(result_matrix) {
  # Get unique IDs
  unique_ids <- unique(result_matrix$id)
  
  # Process each ID
  all_transposed <- list()
  
  for (i in seq_along(unique_ids)) {
    id <- unique_ids[i]
    transposed <- process_single_id(result_matrix, id)
    all_transposed[[i]] <- transposed
  }
  
  # Combine all using rbind (they can have different numbers of columns)
  final_combined <- data.frame()
  
  for (transposed_df in all_transposed) {
    if (nrow(final_combined) == 0) {
      final_combined <- transposed_df
    } else {
      # Align columns by name
      all_cols <- unique(c(names(final_combined), names(transposed_df)))
      
      # Add missing columns to final_combined
      for (col in all_cols) {
        if (!col %in% names(final_combined)) {
          final_combined[[col]] <- NA
        }
      }
      
      # Add missing columns to transposed_df
      for (col in all_cols) {
        if (!col %in% names(transposed_df)) {
          transposed_df[[col]] <- NA
        }
      }
      
      # Reorder columns to match
      transposed_df <- transposed_df[, names(final_combined)]
      
      # Bind rows
      final_combined <- rbind(final_combined, transposed_df)
    }
  }
  
  return(final_combined)
}

transform_and_export <- function(result_matrix, output_file = "results/output_temperature_matrix_t.txt") {
  # Transform to wide format
  final_matrix <- process_all_ids(result_matrix)
  
  # Export to txt file
  write.table(final_matrix,
              file = output_file,
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE,
              na = "NA")
  
  cat("✓ Transformation completed successfully!\n")
  cat("✓ Final matrix dimensions:", nrow(final_matrix), "rows ×", ncol(final_matrix), "columns\n")
  cat("✓ File saved to:", output_file, "\n")
  
  return(final_matrix)
}

transform_and_export(temps_table)
