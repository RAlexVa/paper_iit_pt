rm(list=ls())
library(tidyverse)
library(stringr)


##### Run all of this to read and parse the file #####


chosen_file <- "results/jobs_since_jan8.txt"
{
  # Function to convert SLURM time format to hours
  convert_slurm_time_to_hours <- function(time_str) {
    # Handle NA values
    if (is.na(time_str)) {
      return(NA)
    }
    
    # Replace + with 00
    time_str <- str_replace(time_str, "\\+", "00")
    
    # Check if format includes days (contains a dash)
    if (grepl("-", time_str)) {
      # Format: DD-HH:MM:SS or D-HH:MM:SS
      parts <- str_split(time_str, "-")[[1]]
      days <- as.numeric(parts[1])
      time_part <- parts[2]
      
      # Split the time part
      time_components <- str_split(time_part, ":")[[1]]
      hours <- as.numeric(time_components[1])
      minutes <- as.numeric(time_components[2])
      seconds <- as.numeric(time_components[3])
      
      # Calculate total hours
      total_hours <- days * 24 + hours + minutes / 60 + seconds / 3600
      
    } else {
      # Format: HH:MM:SS
      time_components <- str_split(time_str, ":")[[1]]
      hours <- as.numeric(time_components[1])
      minutes <- as.numeric(time_components[2])
      seconds <- as.numeric(time_components[3])
      
      # Calculate total hours
      total_hours <- hours + minutes / 60 + seconds / 3600
    }
    
    return(total_hours)
  }
  
  
  # Read the file to get all lines
  lines <- readLines(chosen_file)
  
  # The first line is the header, second line has dashes showing column widths
  header_line <- lines[1]
  dash_line <- lines[2]
  
  # Find the positions where dashes start and end (column boundaries)
  # We'll use gregexpr to find sequences of dashes
  dash_positions <- gregexpr("-+", dash_line)[[1]]
  dash_lengths <- attr(dash_positions, "match.length")
  
  # Calculate start and end positions for each column
  col_starts <- dash_positions
  col_ends <- dash_positions + dash_lengths - 1
  
  # Create a specification for read_fwf
  col_positions <- fwf_positions(
    start = col_starts,
    end = col_ends,
    col_names = unlist(strsplit(trimws(header_line), "\\s+"))
  )
  
  # Read the data using fixed width format, skipping the first 2 lines
  jobs_data <- read_fwf(chosen_file, 
                        col_positions = col_positions,
                        skip = 2,
                        trim_ws = TRUE,
                        na = c("", "NA"),
                        col_types = cols(.default = "c"))
  
  #Delete some rows
  jobs_data <- jobs_data |> filter(State %in% c("COMPLETED","TIMEOUT"))
  
  #Delete the decimals from JobID
  jobs_data$JobID <- str_remove(jobs_data$JobID,"[.+].*")
  
  # Now process as before
  jobs_data <- jobs_data %>%
    separate(JobID, into = c("ParentJobID", "ArrayID"), 
             sep = "_", 
             remove = FALSE,
             fill = "right")
  
  #Get the info for MaxRSS
  batch_data <- jobs_data |> 
    filter(JobName=="batch") |> 
    select(JobID,MaxRSS)
  
  #Paste the info into the fila and filter only 1 row per job
  jobs_data <- jobs_data |> select(-MaxRSS) |> 
    left_join(batch_data,by="JobID") |> 
    filter(!(JobName %in% c("batch","extern")))
  
  
  
  # Convert the time in the file to hours
  jobs_data <- jobs_data %>%
    mutate(
      Elapsed_hours = sapply(Elapsed, convert_slurm_time_to_hours),
      Timelimit_hours = sapply(Timelimit, convert_slurm_time_to_hours),
      TotalCPU_hours = sapply(TotalCPU, convert_slurm_time_to_hours),
    )
  
  jobs_data$ParentJobID <- as.numeric(jobs_data$ParentJobID)
  
}
convert_maxrss_to_mb <- function(maxrss_vector) {
  # Function to convert a single MaxRSS value to MB
  sapply(maxrss_vector, function(x) {
    # Handle NA, empty, or actual NA values
    if (is.na(x) || x == "" || x == "NA") {
      return(NA_real_)
    }
    
    # Extract numeric part and suffix
    value <- as.character(x)
    suffix <- substr(value, nchar(value), nchar(value))
    numeric_part <- as.numeric(substr(value, 1, nchar(value) - 1))
    
    # Convert based on suffix
    mb_value <- case_when(
      suffix == "K" ~ numeric_part / 1024,
      suffix == "M" ~ numeric_part,
      suffix == "G" ~ numeric_part * 1024,
      suffix == "T" ~ numeric_part * 1024 * 1024,
      TRUE ~ NA_real_  # For any unexpected format
    )
    
    return(mb_value)
  })
}


#############################################################

#Check start and end
check_start_end <- jobs_data |>
  filter(ParentJobID>6628000) |> 
  select(ParentJobID,Start,End) |> #,ArrayID
  group_by(ParentJobID) |> 
  summarise(min_start=min(Start),max_end=max(End), count=n())
View(check_start_end)


jobs_data$MaxRSS <- convert_maxrss_to_mb(jobs_data$MaxRSS)
#Time taken
summary <- jobs_data |> 
  # filter(ParentJobID %in% chosen_ids) |> 
  filter(ParentJobID>6295000) |>
  group_by(ParentJobID) |> 
  summarise(min_time=min(Elapsed_hours),
            avg_time=mean(Elapsed_hours),
            max_time=max(Elapsed_hours),
            time_limit=max(Timelimit_hours),
            max_RSS_mb=max(MaxRSS,na.rm=T))


state_summary <- jobs_data |> 
  group_by(ParentJobID,State) |> 
  summarise(count=n()) |> ungroup() |> 
  pivot_wider(names_from = State,values_from=count)

summary <- summary |> left_join(state_summary,by="ParentJobID")

# print(summary,n=nrow(summary))
View(summary)



### Checking specific parent id

indiv <- jobs_data |> filter(ParentJobID==6112404)
View(indiv)



