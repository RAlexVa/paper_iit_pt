rm(list=ls())
library(tidyverse)
library(ggplot2)
library(stringr)
library(stargazer)
#### Create plots to include in the PDF
wd_path <- "C:/Users/ralex/Documents/src/paper_iit_pt"
source_path <- "C:/Users/ralex/Documents/src/paper_iit_pt/results"
export_path <- "C:/Users/ralex/Documents/src/paper-adaptive-iit-latex/images/highdim_ex"
export_path_lowdim <- "C:/Users/ralex/Documents/src/paper-adaptive-iit-latex/images/lowdim_ex"
tables_path <- "C:/Users/ralex/Documents/src/paper-adaptive-iit-latex/tables"


##### PLOT FORMAT #####
{
  #Renaming the algorithms
  alg_names <- tibble(alg=c("PT A-IIT(min)","PT A-IIT(sq)","PT-IIT(min)","PT-IIT(sq)"),
                      algorithm=c("MH-mult","A-IIT","RF-MH","IIT"))
  # Colors for each algorithm
  alg_colors <- c(
    "A-IIT" = "#E41A1C",
    "MH-mult" = "#377EB8",
    "IIT" = "#4DAF4A",
    "RF-MH" = "#984EA3"
  )
  
  #Theme for all plots
  my_theme <- theme_minimal() +
    theme(
      # text = element_text(size = 12),           # Base text size
      axis.title = element_text(size = 22),     # Axis labels
      axis.text = element_text(size = 22),      # Axis tick labels
      legend.title = element_text(size = 16),   # Legend title
      legend.text = element_text(size = 14),    # Legend labels
      plot.title = element_text(size = 16)      # Plot title (if you add one)
    )
  theme_set(my_theme)
}
##### FUNCTIONS #####

{
  #To check if the list exists
  check_list <- function(filename, input_ids){
    # Check if file exists
    if (file.exists(filename)) {
      # Read existing file
      existing_list <- readRDS(filename)
      
      # Check if "ids" entry exists in the list
      if ("ids" %in% names(existing_list)) {
        existing_ids <- existing_list$ids
        
        # Compare existing ids with input ids
        # Check if they're identical (same values, same order, same length)
        if (identical(sort(existing_ids), sort(input_ids))) {
          message("IDs are the same. File not modified.");return(TRUE)
        } else {
          message("IDs are different. Update file.");return(FALSE)
        }
      } else {
        # ids entry doesn't exist, save new list
        message("No 'ids' entry found.");return(FALSE);
      }
    } else {
      # File doesn't exist, create it
      message("File doesn't exist. Create new file.");return(FALSE);
    }
  }
  #To create the list with the consolidated information
  create_list <- function(chosen_ids,chosen_dim){
    # chosen_dim <- "highdim"
    consolidated_list <- list()
    if(chosen_dim=="highdim"){
      ##### Define IDs and dimension #####
      
      parameters <- read_csv(file.path(wd_path,"inputs/simulation_details_highd.csv"), col_types = cols())
      
      #List of files in the results folder
      data_sum <- tibble(file_names=list.files(path = source_path, pattern = "^sim_.*\\.Rds")) |> 
        mutate(id=as.numeric(str_extract(file_names, "id_([0-9]+)",group=1)),
               sim_id=str_extract(file_names, "id_[0-9]+_([0-9]+)\\.Rds",group=1),
               dim=str_extract(file_names, "(?<=sim_)[^_]+(?=_id)")) |> 
        filter(dim==chosen_dim) |> 
        left_join(parameters, by="id")|> 
        filter(id %in% chosen_ids)
      
      
      ##### Create dataset to store results#####
      
      #Check which columns are full of NAs to exclude them
      na_temps <- data_sum |> 
        select(matches("^t\\d+$")) |>  # Select columns of the form tX
        summarize(across(everything(), ~ all(is.na(.x)))) |>   # Check if all values are NA
        pivot_longer(everything(), names_to = "column", values_to = "is_all_na") |>   # Convert to long format
        filter(is_all_na) |>   # Keep only columns where is_all_na is TRUE
        pull(column) 
      
      #Delete the columns for temperatures with NA
      data_sum <- data_sum |> select(-all_of(na_temps))
      
      temporal_dataset <- data_sum |> 
        select(id,matches("^t\\d+$")) |> 
        group_by(id) |> 
        slice(1) |> ungroup()
      #Compute max. number of replicas in the chosen IDs
      max_replicas <- ncol(temporal_dataset)-1
      
      cat("Report on # of temperatures")
      print(
        temporal_dataset |> rowwise() |> 
          mutate(temp_count=max_replicas-sum(is.na(c_across(-id)))) |> 
          select(id,temp_count) 
      )
      
      ##### Define datasets to store results#####
      distances <- tibble(alg=character(),sim=numeric(),temperature=numeric(),iteration=numeric(),mode=character())
      round_trip <- tibble(alg=character(),sim=numeric(),replica=numeric(),round_trips=numeric());
      swap_rate <- tibble(alg=character(),sim=numeric(),replica=numeric(),swap_rate=numeric());
      iterations <- tibble(alg=character(),sim=numeric(),replica=numeric(),iterations=numeric());
      log_bounds <- tibble(alg=character(),sim=numeric(),replica=numeric(),log_bound=numeric());
      rf_replicas <- tibble(id=numeric(),rf_reps=numeric());
      temperatures_matrix <- tibble(id=numeric(),temp_id=numeric(),temperature=numeric());
      
      time_taken <- as.data.frame(matrix(nrow=0,ncol=2));colnames(time_taken) <- c("alg","time")
      
      ##### FOR loop to read .Rds files and get results#####
      # Start creating datasets with information
      for(i in 1:nrow(data_sum)){
        data <- readRDS(file.path(source_path,data_sum[i,1]))
        selected_id <- data_sum |> slice(i)|> pull(id)#Extract ID
        tot_sim <- data_sum |> slice(i)|> pull(simulations) #Should just be 1
        algorithm <- data_sum |> slice(i) |> pull(algorithm)
        tot_iter <- data_sum |> slice(i) |> pull(iterations)
        tot_swap <- data_sum |> slice(i) |> pull(total_swap)
        interswap <- data_sum |> slice(i) |> pull(interswap)
        temperatures <- as.numeric(data_sum |> slice(i) |> select(matches("^t\\d{1,2}$")))
        temperatures <- temperatures[!is.na(temperatures)]# all the temperatures are in order and consecutive in the CSV file

        num_modes <- data_sum |> slice(i) |> pull(num_modes)
        #Add the BF to the algorithm name
        if(algorithm=="PT_A_IIT"){algorithm <- paste0("PT A-IIT(",data_sum |> slice(i)|> pull(bf),")")}
        if(algorithm=="PT_IIT_Z"){algorithm <- paste0("PT-IIT(",data_sum |> slice(i)|> pull(bf),")")}
        
        #Temperatures_matrix
        temperatures_temporal <- tibble(id=selected_id,temp_id=1:length(temperatures),temperature=temperatures,alg=algorithm);
        temperatures_matrix <- rbind(temperatures_matrix,temperatures_temporal)
        
        
        ###The simulation ID is in the name
        sim_id <- as.numeric(data_sum |> slice(i) |> pull(sim_id))
        
        
        ##### Optional add the ID of the simulation into the name of the algorithm
        print(paste0(selected_id," ",sim_id," ",algorithm))
        algorithm <- paste0(algorithm,"(",selected_id,")")
        
        #For each selected_id extract the number of Rejection Free replicas
        temp <- data_sum |> slice(i) |> pull(temps_rf)
        if(!(selected_id %in% rf_replicas$id)){
          rf_replicas <- rf_replicas |> add_row(id=selected_id,rf_reps=temp)
        }
        
        
        #Extract time
        temp <- as.data.frame(data[["time_taken"]])
        if(!is_empty(temp)){
          temp$alg <- algorithm
          temp <- temp |> select(alg,everything())
          colnames(temp) <- c("alg","time")
          time_taken <- rbind(time_taken,temp)
        }
        
        if(chosen_dim=="highdim"){ ## Extract distance to modes and time to reach it
          ### Specific extractions for highdim example
          p <- data_sum |> slice(i) |> pull(p);
          
          ##### Extract distance to modes
          output_name <- paste0("distance_modes")
          output_time <- paste0("time_modes")
          
          ### Extract minimum distances  
          temp_m <- as.data.frame(data[[output_name]]) 
          colnames(temp_m) <- round(temperatures,2) #column name is the temperature value
          #Create columns to identify
          temp_m$alg <- algorithm #Column with name of the algorithm
          temp_m$mode <- paste0("m",1:num_modes) #Column indicating the mode
          temp_m$sim <- sim_id #Column indicating the number of simulation
          
          #Pivot longer
          temp_m <- temp_m |> pivot_longer(-(alg:sim),names_to="temperature",values_to = "min_dist")
          
          ### Extract time to reach minimum distances
          temporal_time <- as.data.frame(data[[output_time]])
          colnames(temporal_time) <- round(temperatures,2)
          #Create columns to identify
          temporal_time$alg <- algorithm
          temporal_time$mode <- paste0("m",1:num_modes)
          temporal_time$sim <- sim_id
          #Pivot longer
          temporal_time <- temporal_time |> pivot_longer(-(alg:sim),names_to="temperature",values_to = "time_find")
          #Join distances and times
          temp_join <- left_join(temp_m,temporal_time,by=c("alg","mode","sim","temperature"))
          #Add to the report
          distances <- rbind(distances,temp_join)
          
        }
        
        
        if(!is_empty(data[["round_trips"]]) && ncol(data[["round_trips"]])>1){
          ### Extract number of round trips rate
          temp <- as.data.frame(data[["round_trips"]])
          colnames(temp) <- 1:ncol(temp)
          temp$sim <- sim_id ### IMPORTANT: This works because there's only 1 simulation in the output
          temp$alg <- algorithm
          #Pivot longer so we can have different number of replicas
          temp <- temp |> pivot_longer(-(alg:sim),names_to="replica",values_to = "round_trips")
          round_trip <- rbind(round_trip,temp)
          
          ### Extract replica swap rate
          temp <- as.data.frame(data[["swap_rate"]])
          colnames(temp) <- 1:ncol(temp)
          temp$sim <- sim_id
          temp$alg <- algorithm
          temp <- temp |> pivot_longer(-(alg:sim),names_to="replica",values_to = "swap_rate")
          swap_rate <- rbind(swap_rate,temp)
        }
        if(!is_empty(data[["total_iter"]])&& ncol(data[["total_iter"]])>1){ 
          ### Extract avg. number of iterations interswap
          #Considering that we only have 1 simulation per read dataset
          temp_single <- data[["total_iter"]][,,1]
          final_swap <- data[["final_swap"]]
          if(!is_empty(final_swap)){
            cutoff <- final_swap
          }else{
            check_rows <- rowSums(temp_single)
            cutoff <- min(which(check_rows==0))-1#This could be INF
            cutoff <- min(cutoff,length(check_rows))#In case we don't have any zeroes = didn't break the for loop
          }
          #Delete the rows without info
          temp_single <- temp_single[1:cutoff,]
          temp <- data.frame(matrix(colSums(temp_single)/nrow(temp_single),nrow=1))
          colnames(temp) <- 1:length(temp)
          temp$sim <- sim_id
          temp$alg <- algorithm
          temp <- temp |> pivot_longer(-(alg:sim),names_to="replica",values_to = "iterations")
          iterations <- rbind(iterations,temp)
        }
        ### Extract log-bounds
        if(!is_empty(data[["final_bounds"]])){
          temp <- as.data.frame(data[["final_bounds"]])# IMPORTANT: Here we assume there's only 1 column
          colnames(temp) <- "log_bound"
          temp$replica <- 1:length(temperatures)
          temp$sim <- sim_id
          temp$alg <- algorithm
          temp <- temp |> select(alg,sim,replica,log_bound)
          log_bounds <- rbind(log_bounds,temp)
        }
        
      }
      
      # Create list with all the data
      
      consolidated_list[["ids"]] <- chosen_ids;
      consolidated_list[["distances"]] <- distances;
      consolidated_list[["round_trip"]] <- round_trip;
      consolidated_list[["swap_rate"]] <- swap_rate;
      consolidated_list[["iterations"]] <- iterations;
      consolidated_list[["log_bounds"]] <- log_bounds;
      consolidated_list[["rf_replicas"]] <- rf_replicas;
      consolidated_list[["temperatures_matrix"]] <- temperatures_matrix;
      consolidated_list[["time_taken"]] <- time_taken;
    }
    
    if(chosen_dim=="lowdim"){
      parameters <- read_csv(paste0("inputs/simulation_details_lowd.csv"), col_types = cols())
      
      #List of files in the results folder
      data_sum <- tibble(file_names=list.files(path = "results", pattern = "^sim_.*\\.Rds")) |> 
        mutate(id=as.numeric(str_extract(file_names, "id_([0-9]+)",group=1)),
               sim_id=str_extract(file_names, "id_[0-9]+_([0-9]+)\\.Rds",group=1),
               dim=str_extract(file_names, "(?<=sim_)[^_]+(?=_id)")) |> 
        filter(dim==chosen_dim) |> 
        left_join(parameters, by="id")|> 
        filter(id %in% chosen_ids)
      
      ##### Create dataset to store results#####
      
      #Check which columns are full of NAs to exclude them
      na_temps <- data_sum |> 
        select(matches("^t\\d+$")) |>  # Select columns of the form tX
        summarize(across(everything(), ~ all(is.na(.x)))) |>   # Check if all values are NA
        pivot_longer(everything(), names_to = "column", values_to = "is_all_na") |>   # Convert to long format
        filter(is_all_na) |>   # Keep only columns where is_all_na is TRUE
        pull(column) 
      
      #Delete the columns for temperatures with NA
      data_sum <- data_sum |> select(-all_of(na_temps))
      
      
      temporal_dataset <- data_sum |> 
        select(id,matches("^t\\d+$")) |> 
        group_by(id) |> 
        slice(1) |> ungroup()
      #Compute max. number of replicas in the chosen IDs
      max_replicas <- ncol(temporal_dataset)-1
      
      cat("Report on # of temperatures")
      print(
        temporal_dataset |> rowwise() |> 
          mutate(temp_count=max_replicas-sum(is.na(c_across(-id)))) |> 
          select(id,temp_count) 
      )
      
      ##### Define datasets to store results#####
      round_trip <- tibble(alg=character(),sim=numeric(),replica=numeric(),round_trips=numeric());
      swap_rate <- tibble(alg=character(),sim=numeric(),replica=numeric(),swap_rate=numeric());
      iterations <- tibble(alg=character(),sim=numeric(),replica=numeric(),iterations=numeric());
      log_bounds <- tibble(alg=character(),sim=numeric(),replica=numeric(),log_bound=numeric());
      temperatures_matrix <- tibble(id=numeric(),temp_id=numeric(),temperature=numeric());
      time_visit <- tibble(alg=character(),sim=numeric(),mode=character(),time=numeric())
      tvd_report <- tibble(alg=character(),sim=numeric(),measurement=numeric(),tvd=numeric(),time=numeric());
      prob_modes <- tibble(alg=character(),sim=numeric(),mode=character(),prob=numeric());
      rf_replicas <- tibble(id=numeric(),rf_reps=numeric());
      
      mod1 = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1);
      mod2 = c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0);
      mod3 = c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1);
      mod4 = c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0);
      mod5 = c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1);
      mod6 = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1);
      mod7 = c(0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0);
      
      matrix_mode <- cbind(mod1,mod2,mod3,mod4,mod5,mod6,mod7)
      index_mode <- c(65536,21846,43691,256,65281,32770,385)
      
      time_taken <- as.data.frame(matrix(nrow=0,ncol=2));colnames(time_taken) <- c("alg","time")
      
      for(i in 1:nrow(data_sum)){
        data <- readRDS(file.path("results",data_sum[i,1]))
        selected_id <- data_sum |> slice(i)|> pull(id)#Extract ID
        tot_sim <- data_sum |> slice(i)|> pull(simulations) #Should just be 1
        algorithm <- data_sum |> slice(i) |> pull(algorithm)
        tot_iter <- data_sum |> slice(i) |> pull(iterations)
        tot_swap <- data_sum |> slice(i) |> pull(total_swap)
        interswap <- data_sum |> slice(i) |> pull(interswap)
        temperatures <- as.numeric(data_sum |> slice(i) |> select(matches("^t\\d{1,2}$")))
        temperatures <- temperatures[!is.na(temperatures)]# all the temperatures are in order and consecutive in the CSV file
        #If the first temperature is negative then we can have a different definition of temperatures
        if(sign(temperatures[1])==-1){
          temp_from <- temperatures[2]
          temp_to <- temperatures[3]
          temp_total_number <- temperatures[4]
          
          temperatures <- seq(temp_from,temp_to,length.out=temp_total_number)
          
        }
        
       
          check_model <- data_sum |> slice(i) |> pull(model)
          if(check_model=="bimodal"){num_modes <- 2; spec_matrix_mode <- matrix_mode[,2:3];spec_index_mode <- index_mode[2:3];}
          if(check_model=="7_mode"){num_modes <- 7;spec_matrix_mode <- matrix_mode;spec_index_mode <- index_mode;}
        
        
        #Add the BF to the algorithm name
        if(algorithm=="PT_A_IIT"){algorithm <- paste0("PT A-IIT(",data_sum |> slice(i)|> pull(bf),")")}
        if(algorithm=="PT_IIT_Z"){algorithm <- paste0("PT-IIT(",data_sum |> slice(i)|> pull(bf),")")}
        
          #Temperatures_matrix
          temperatures_temporal <- tibble(id=selected_id,temp_id=1:length(temperatures),temperature=temperatures,alg=algorithm);
          temperatures_matrix <- rbind(temperatures_matrix,temperatures_temporal)
          
          
          
        ###The simulation ID is in the name
        sim_id <- as.numeric(data_sum |> slice(i) |> pull(sim_id))
        
        
        ##### Optional add the ID of the simulation into the name of the algorithm
        print(paste0(selected_id," ",sim_id," ",algorithm))
        algorithm <- paste0(algorithm,"(",selected_id,")")
        
        
        
        
        #Extract time
        temp <- as.data.frame(data[["time_taken"]])
        if(!is_empty(temp)){
          temp$alg <- algorithm
          temp <- temp |> select(alg,everything())
          colnames(temp) <- c("alg","time")
          time_taken <- rbind(time_taken,temp)
        }
        
        
        if(chosen_dim=="lowdim"){
          
          #Extract time to visit the modes
          temp <- as.data.frame(data[["time_visit"]])
          temp <- temp |> slice(spec_index_mode)
          colnames(temp) <- "time"
          temp$sim <- sim_id ### IMPORTANT: This works because there's only 1 simulation in the output
          temp$alg <- algorithm
          temp$mode <- paste0("m",1:num_modes) #Column indicating the mode
          temp$sim <- sim_id #Column indicating the number of simulation
          
          time_visit <- rbind(time_visit,temp)
          
          #Extract evolution of TVD
          temp <- as.data.frame(data[["tvd_report"]])
          colnames(temp) <- "tvd"
          temp$sim <- sim_id ### IMPORTANT: This works because there's only 1 simulation in the output
          temp$alg <- algorithm
          temp$measurement <- 1:nrow(temp)
          
          
          temp_time <- as.numeric(data[["tvd_time_report"]])
          temp$time <- temp_time
          temp <- temp |> filter(tvd>0)
          
          tvd_report <- rbind(tvd_report,temp)
          
          #Extract estimated probability of the modes    
          temp <- as.data.frame(data[["est_pi"]])
          
          temp <- temp |> slice(spec_index_mode)
          colnames(temp) <- "prob"
          temp$sim <- sim_id ### IMPORTANT: This works because there's only 1 simulation in the output
          temp$alg <- algorithm
          temp$mode <- paste0("m",1:num_modes) #Column indicating the mode
          temp$sim <- sim_id #Column indicating the number of simulation
          
          prob_modes <- rbind(prob_modes,temp)
          
        }
        
        #For each selected_id extract the number of Rejection Free replicas
        temp <- data_sum |> slice(i) |> pull(temps_rf)
        if(!(selected_id %in% rf_replicas$id)){
          rf_replicas <- rf_replicas |> add_row(id=selected_id,rf_reps=temp)
        }
        
        if(!is_empty(data[["round_trips"]]) && ncol(data[["round_trips"]])>1){
          ### Extract number of round trips rate
          temp <- as.data.frame(data[["round_trips"]])
          colnames(temp) <- 1:ncol(temp)
          temp$sim <- sim_id ### IMPORTANT: This works because there's only 1 simulation in the output
          temp$alg <- algorithm
          #Pivot longer so we can have different number of replicas
          temp <- temp |> pivot_longer(-(alg:sim),names_to="replica",values_to = "round_trips")
          round_trip <- rbind(round_trip,temp)
          
          ### Extract replica swap rate
          temp <- as.data.frame(data[["swap_rate"]])
          colnames(temp) <- 1:ncol(temp)
          temp$sim <- sim_id
          temp$alg <- algorithm
          temp <- temp |> pivot_longer(-(alg:sim),names_to="replica",values_to = "swap_rate")
          swap_rate <- rbind(swap_rate,temp)
        }
        if(!is_empty(data[["total_iter"]])&& ncol(data[["total_iter"]])>1){
          ### Extract avg. number of iterations interswap
          #Considering that we only have 1 simulation per read dataset
          temp_single <- data[["total_iter"]][,,1]
          final_swap <- data[["final_swap"]]
          if(!is_empty(final_swap)){
            cutoff <- final_swap
          }else{
            check_rows <- rowSums(temp_single)
            cutoff <- min(which(check_rows==0))-1#This could be INF
            cutoff <- min(cutoff,length(check_rows))#In case we don't have any zeroes = didn't break the for loop
          }
          #Delete the rows without info
          temp_single <- temp_single[1:cutoff,]
          temp <- data.frame(matrix(colSums(temp_single)/nrow(temp_single),nrow=1))
          colnames(temp) <- 1:length(temp)
          temp$sim <- sim_id
          temp$alg <- algorithm
          temp <- temp |> pivot_longer(-(alg:sim),names_to="replica",values_to = "iterations")
          iterations <- rbind(iterations,temp)
        }
        ### Extract log-bounds
        if(!is_empty(data[["final_bounds"]])){
          temp <- as.data.frame(data[["final_bounds"]])# IMPORTANT: Here we assume there's only 1 column
          colnames(temp) <- "log_bound"
          temp$replica <- 1:length(temperatures)
          temp$sim <- sim_id
          temp$alg <- algorithm
          temp <- temp |> select(alg,sim,replica,log_bound)
          log_bounds <- rbind(log_bounds,temp)
        }
        
      }# End of for loop to read .Rds files
      
      # Create list with all the data
      
      consolidated_list[["ids"]] <- chosen_ids;
      consolidated_list[["round_trip"]] <- round_trip;
      consolidated_list[["swap_rate"]] <- swap_rate;
      consolidated_list[["iterations"]] <- iterations;
      consolidated_list[["log_bounds"]] <- log_bounds;
      consolidated_list[["rf_replicas"]] <- rf_replicas;
      consolidated_list[["temperatures_matrix"]] <- temperatures_matrix;
      consolidated_list[["time_taken"]] <- time_taken;
      consolidated_list[["tvd_report"]] <-tvd_report;
      consolidated_list[["prob_modes"]] <-prob_modes;
      consolidated_list[["time_visit"]] <-time_visit;
      
    }# End of IF statement for lowdim
    
    return(consolidated_list)
    
  }
  #Use the two functions above to import the list with information.
  create_plot_input <- function(filename,input_ids,chosen_dim){
    file_to_create <- file.path(source_path,paste0(filename,".Rds"))
    if(check_list(file_to_create,input_ids)){#If the list already exists and matches
      lista_creada <- readRDS(file_to_create);
    }else{#If the list doesn't exist or doesn't match
      lista_creada <- create_list(input_ids,chosen_dim)
      saveRDS(lista_creada,file_to_create);
    }
    return(lista_creada);
  }
  
  speed_plot_last <- function(the_list, subset_ids){
    distances <- the_list[["distances"]]
    if(!missing(subset_ids)){
      #Subset the distances dataset only for the specified ids
      distances$id <- str_extract(distances$alg, "\\(\\d+\\)") |> 
        str_remove_all("[()]")
      distances <- distances |> 
        filter(id %in% subset_ids) |> 
        select(-id)
    }
    
    
    dist_t1 <- distances |> mutate(temperature=as.numeric(temperature)) |> 
      group_by(alg) |> 
      summarise(max_temp=max(temperature)) |> 
      right_join(distances,by="alg") |> 
      filter(temperature==max_temp) |> 
      select(-max_temp) |> 
      ungroup()
    
    dist_t1_times <- dist_t1  |> 
      mutate(time_fixed=ifelse(min_dist>0,Inf,time_find)) |> 
      select(-min_dist,-time_find) |> 
      pivot_longer(cols=c(time_fixed),names_to="variable",values_to="measure") |> 
      pivot_wider(names_from=mode,values_from=measure) |> 
      rowwise() |> 
      mutate(first_time = min(across(matches("^m\\d+$"))),
             last_time = max(across(matches("^m\\d+$"))))
    dist_t1_times[dist_t1_times==Inf] <- NA ## replace Infinities for NAs
    {
      report_last <- dist_t1_times |> 
        select(alg,last_time)|> 
        rename(time_visit=last_time)
      # fit <- survfit(Surv(time_visit,rep(1,nrow(report_last)))~alg,data=report_last)
      
      
      # (plot_surv_mode <- ggsurvplot(fit,
      #                               data=report_last,
      #                               fun="event",
      #                               palette = "Set1",    # Color palette
      #                               xlab = "Time (seconds)",
      #                               ylab = "Prop. of simulations visiting all modes",
      #                               legend.title = "Algorithm",
      #                               # break.time.by = time_br,
      #                               font.x = 25,        # X-axis label font size
      #                               font.y = 25,        # Y-axis label font size
      #                               font.tickslab = 18, # Axis tick labels (numbers) font size
      #                               font.legend = 20,
      #                               size=2,
      #                               conf.int = FALSE,
      #                               censor = TRUE))   # Legend text font size)
      
    }
    #Modify names of algorithm
    report_last$alg <- str_remove(report_last$alg,"\\([^)]*\\)$")
    report_last <- report_last |> left_join(alg_names,by="alg") |> ungroup()
    
    # plot_surv_mode <- ggplot(report_last,aes(x=time_visit,color=algorithm))+
    #   stat_ecdf(linewidth=1.6)+
    #   scale_color_manual(values = alg_colors, name = "Algorithm") +
    #   scale_y_continuous(labels = scales::percent_format()) +
    #   labs(x = "seconds", y = "%")
    
    
    # Calculate eCDF manually for each algorithm
    ecdf_data <- report_last |> 
      filter(!is.na(time_visit)) |> 
      group_by(alg) |> 
      arrange(time_visit) |> 
      mutate(
        ecdf_value = (row_number()) / n()
      ) |> 
      ungroup()
    
    # Plot with geom_step
    plot_surv_mode <- ggplot(ecdf_data, aes(x = time_visit, y = ecdf_value, color = algorithm)) +
      geom_step(linewidth=1.6) +
      scale_color_manual(values = alg_colors, name = "Algorithm") +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(x = "seconds", y = "% of simulations")
    
    
    
    
    
    
    return(plot_surv_mode)
  }
  
  violin_plot_last <- function(the_list,add_info,subset_ids,subset_rfs, return_dataset=F){
    distances <- the_list[["distances"]]
    distances$id <- as.numeric(str_extract(distances$alg, "\\(\\d+\\)") |> 
                                 str_remove_all("[()]"))
    
    distances <- distances |> left_join(add_info,by="id")
    
    #Subset the distances dataset only for the specified ids
    if(!missing(subset_ids)){
      distances <- distances |> 
        filter(id %in% subset_ids) 
    }
    if(!missing(subset_rfs)){
      #Subset the distances dataset only for the specified ids
      distances <- distances |> 
        filter(rf %in% subset_rfs) 
    }
    distances <- distances |> select(-id);
    
    dist_t1 <- distances |> mutate(temperature=as.numeric(temperature)) |> 
      group_by(alg) |> 
      summarise(max_temp=max(temperature)) |> 
      right_join(distances,by="alg") |> 
      filter(temperature==max_temp) |> 
      select(-max_temp) |> 
      ungroup()
    
    dist_t1_times <- dist_t1  |> 
      mutate(time_fixed=ifelse(min_dist>0,Inf,time_find)) |> 
      select(-min_dist,-time_find) |> 
      pivot_longer(cols=c(time_fixed),names_to="variable",values_to="measure") |> 
      pivot_wider(names_from=mode,values_from=measure) |> 
      rowwise() |> 
      mutate(first_time = min(across(matches("^m\\d+$"))),
             last_time = max(across(matches("^m\\d+$"))))
    dist_t1_times[dist_t1_times==Inf] <- NA ## replace Infinities for NAs
    
    report_last <- dist_t1_times |> 
      select(alg,rf,last_time)|> 
      rename(time_visit=last_time)
    
    #Modify names of algorithm
    report_last$alg <- str_remove(report_last$alg,"\\([^)]*\\)$")
    report_last <- report_last |> left_join(alg_names,by="alg") |> ungroup()
    
    plot_surv_mode <- ggplot(report_last,aes(x=factor(rf),fill=algorithm, y=time_visit))+
      geom_violin()+
      scale_fill_manual(values = alg_colors, name = "Algorithm") +
      labs(x = "Rejection Free replicas", y = "seconds")
    
    ret_list=list();
    
    if(return_dataset){
      return(report_last)
    }else{
      return(plot_surv_mode) 
    }
    
  }
  
  export_plot <- function(the_plot,plot_name,h_l_dim,dimension){
    if(h_l_dim !="lowdim" && h_l_dim != "highdim"){stop("Incorrect dimension. Should be lowdim or highdim")}
    
    if(h_l_dim=="lowdim"){
      if(missing(dimension)){
        plot_name <- paste0(plot_name,"_ld.jpg")
        file_name <- file.path(export_path_lowdim,plot_name)
        
      }else{
        plot_name <- paste0(plot_name,"_ld_",dimension,".jpg")
        file_name <- file.path(export_path_lowdim,plot_name)
      }

    }
    
    if(h_l_dim=="highdim"){
      if(missing(dimension)){
        plot_name <- paste0(plot_name,"_hd.jpg")
        file_name <- file.path(export_path,plot_name)
        
      }else{
        plot_name <- paste0(plot_name,"_hd_",dimension,"k.jpg")
        file_name <- file.path(export_path,plot_name)
      }
      
    }
    
    

    
    jpeg(file_name,width=1200,height = 600,pointsize = 30)
    print(the_plot)
    dev.off()
    message("Successfully created plot: ",plot_name)
  }
  
  # To create swap rate table in latex
  create_swap_rate_table <- function(lll,h_l_dim,dimension){
    if(h_l_dim !="lowdim" && h_l_dim != "highdim"){stop("Incorrect dimension. Should be lowdim or highdim")}
    
    
    swap_rate <- lll$swap_rate
    swap_rate$replica <- as.numeric(swap_rate$replica)
    
    swap_rate <- fix_alg_name(swap_rate)
    
    sr_report <- swap_rate |>
      group_by(algorithm,replica) |> 
      summarise(avg.sr=mean(swap_rate,na.rm=T)) |> ungroup() |> 
      arrange(algorithm,replica)
    
    #Horizontal
    # horizontal_report_sr <- sr_report |> 
    #   pivot_wider(names_from = replica,values_from = avg.sr)
    
    #Vertical
    vertical_report_sr <- sr_report  |> 
      pivot_wider(names_from = algorithm,values_from = avg.sr)
    
    # boxplot_sr <- swap_rate |> 
    #   ggplot(aes(x=factor(replica),y=swap_rate))+
    #   geom_boxplot()+
    #   geom_hline(yintercept=0.234,col="red")+
    #   facet_wrap(~algorithm)
    
    # table_sr_report <- horizontal_report_sr |>
    #   mutate(across(-algorithm,round,digits=4))
    table_sr_report <- vertical_report_sr |>
      mutate(across(-replica,round,digits=4))
    
    #Export table with averages
    # stargazer(table_sr_report,
    #           summary=F,
    #           type="latex",
    #           rownames=F,
    #           title=paste0("Average swap rate between replica $\\beta_{i}$ and  $\\beta_{i+1}$"),
    #           label=paste0("tab:swap_rate_",dimension,"k"),
    #           out=file.path(tables_path,paste0("avg_swap_rate_dim",dimension,".tex")))
    
    #Particular changes depending on model
    table_file_name <- ""
    table_caption <- ""
    table <- label <- ""
    if(h_l_dim=="lowdim"){
      table_file_name <- paste0("avg_swap_rate_ld_", dimension, ".tex");
      if(dimension=="7_mode"){dimension_patch <-"7 modes" }
      if(dimension=="bimodal"){dimension_patch <-"bimodal" }
      table_caption <- paste0("Average swap rate between replica $\\beta_{i}$ and $\\beta_{i+1}$ for the ",dimension_patch," model")
      table_label <- paste0("tab:swap_rate_", dimension)
      }
    if(h_l_dim=="highdim"){
      table_file_name <- paste0("avg_swap_rate_hd_dim", dimension, ".tex");
      table_caption <- paste0("Average swap rate between replica $\\beta_{i}$ and $\\beta_{i+1}$ in space of dimension ",dimension,"000")
      table_label <- paste0("tab:swap_rate_", dimension, "k")
      }
    
    #Export modifying the name of rows
    # Generate stargazer output as a character vector
    latex_code <- stargazer(table_sr_report,
                            summary = FALSE,
                            type = "latex",
                            rownames = FALSE,
                            title = table_caption,
                            label = table_label)
    
    # Replace numbers in replica with beta and subindex
    for (i in table_sr_report$replica) {
      latex_code <- gsub(paste0("^", i, " &"), 
                         paste0("$\\\\beta_{", i, "}$ &"), 
                         latex_code)
    }
    
    
    
    # Write to file
    writeLines(latex_code, 
               file.path(tables_path, table_file_name))
    
    
  }
  
  
  # To create latex table of average number of iterations per replica
  create_table_iterations <- function(the_list,h_l_dim,dimension,subset_ids){
    dataset <- the_list[["iterations"]]
    
    dataset$id <- as.numeric(str_extract(dataset$alg, "\\(\\d+\\)") |> 
                               str_remove_all("[()]")) #Create ID variable
    dataset$replica <- as.numeric(dataset$replica)
    dataset <- fix_alg_name(dataset)
    
    if(!missing(subset_ids)){dataset <- dataset |> filter(id %in% subset_ids)}
    
    summary_report <- dataset |> 
      group_by(algorithm,id,replica) |> 
      summarise(avg.iter=mean(iterations)) |> 
      ungroup()
    
    summary_report <- summary_report |> 
      select(algorithm,replica,avg.iter) |> 
      pivot_wider(names_from=algorithm,values_from=avg.iter) |> 
      mutate(across(-replica, ~ round(.x, 3)))
    
    
    #Particular changes depending on model
    table_file_name <- ""
    table_label <- ""
    table_caption <- paste0("Average number of iterations between replica swaps")
    if(h_l_dim=="lowdim"){
      table_file_name <- paste0("avg_iter_ld_", dimension, ".tex");
      
      table_label <-paste0("tab:avg_iter_", dimension)
    }
    if(h_l_dim=="highdim"){
      table_file_name <- paste0("avg_iter_hd_dim", dimension, ".tex");
      table_label <-paste0("tab:avg_iter_", dimension, "k")
    }
    
    #Export modifying the name of rows
    # Generate stargazer output as a character vector
    latex_code <- stargazer(summary_report,
                            summary = FALSE,
                            type = "latex",
                            rownames = FALSE,
                            title = table_caption,
                            label = table_label)
    
    # Replace numbers in replica with beta and subindex
    for (i in summary_report$replica) {
      latex_code <- gsub(paste0("^", i, " &"), 
                         paste0("$\\\\beta_{", i, "}$ &"), 
                         latex_code)
    }
    
    
    
    # Write to file
    writeLines(latex_code, 
               file.path(tables_path, table_file_name))
    
    
    
  }
  
  # To create latex table of the values of the betas used in the replicas
  create_table_temperatures <- function(the_list,h_l_dim,dimension,subset_ids){
    dataset <- the_list[["temperatures_matrix"]]
    
    
    dataset <- fix_alg_name(dataset)
    
    if(!missing(subset_ids)){dataset <- dataset |> filter(id %in% subset_ids)}
    
    dataset <- unique(dataset)
    
    table_temperatures <- dataset |> 
      select(algorithm,temperature,temp_id) |> 
      pivot_wider(names_from = algorithm,values_from=temperature)
    
    
    #Particular changes depending on model
    table_file_name <- ""
    table_label <- ""
    
    if(h_l_dim=="lowdim"){
      table_file_name <- paste0("temperatures_ld_", dimension, ".tex");
      table_label <-paste0("tab:temperatures_", dimension)
      if(dimension=="7_mode"){dimension_patch <-"7 modes" }
      if(dimension=="bimodal"){dimension_patch <-"bimodal" }
      table_caption <- paste0("Inverse temperature used in the ",dimension_patch," model")
    }
    if(h_l_dim=="highdim"){
      table_file_name <- paste0("temperatures_hd_dim", dimension, ".tex");
      table_label <-paste0("tab:temperatures_", dimension, "k")
      table_caption <- paste0("Inverse temperature used in the problem of dimension ",dimension,"k")
    }
    
    #Export modifying the name of rows
    # Generate stargazer output as a character vector
    latex_code <- stargazer(table_temperatures,
                            summary = FALSE,
                            type = "latex",
                            rownames = FALSE,
                            title = table_caption,
                            label = table_label)
    
    # Replace numbers in replica with beta and subindex
    for (i in table_temperatures$temp_id) {
      latex_code <- gsub(paste0("^", i, " &"), 
                         paste0("$\\\\beta_{", i, "}$ &"), 
                         latex_code)
    }
    
    
    
    # Write to file
    writeLines(latex_code, 
               file.path(tables_path, table_file_name))
  }
  
  
### Specifically for lowdim  
  # To create plot of Total Variation Distance
  plot_tvd <- function(the_list,filter_measurement=50,filter_time=1500, subset_ids){
    tvd_report <- the_list[["tvd_report"]] #Extract tvd_report
    tvd_report$id <- as.numeric(str_extract(tvd_report$alg, "\\(\\d+\\)") |> 
                                  str_remove_all("[()]")) #Create ID variable
    
    tvd_report <- fix_alg_name(tvd_report)
    tvd_summary <- tvd_report |> 
      group_by(id,algorithm,measurement) |> 
      summarise(avg_tvd=mean(tvd),avg_time=mean(time)) |> 
      ungroup()
    
    if(!missing(subset_ids)){#If we want to subset IDs
      tvd_summary <- tvd_summary |> filter(id %in% subset_ids)
    }
    
    tvd_plot <- tvd_summary |> 
      filter(measurement>filter_measurement) |>
      filter(avg_time<filter_time) |> 
      ggplot(aes(x=avg_time,y=avg_tvd,color=algorithm))+
      geom_line(linewidth=1.2, alpha=0.4)+
      geom_point(aes(x=avg_time,y=avg_tvd),size=0.4,alpha=0.5)+
      scale_color_manual(values = alg_colors, name = "Algorithm")+
      labs(x = "seconds", y = "TVD")
    
    return(tvd_plot)
    
  }
  
  # To create plot of speed to mode
  speed_mode_lowdim <- function(the_list,subset_ids){
    time_visit <- the_list[["time_visit"]]
    time_visit$id <- as.numeric(str_extract(time_visit$alg, "\\(\\d+\\)") |> 
                                  str_remove_all("[()]"))
    time_visit <- fix_alg_name(time_visit)
    
    if(!missing(subset_ids)){
      time_visit <- time_visit |> filter(id %in% subset_ids)
    }
    
    time_visit_report <- time_visit |> 
      pivot_wider(names_from=mode, values_from=time) |> 
      mutate(count_modes=rowSums(!is.na(across(starts_with("m"))))) |> 
      rowwise() |> 
      mutate(first_time = min(across(matches("^m\\d+$")),na.rm=T),
             last_time = max(across(matches("^m\\d+$")),na.rm=T))
    #Manually calculate ECDF
    ecdf_data <- time_visit_report |> 
      select(algorithm,last_time) |> 
      group_by(algorithm) |> 
      arrange(last_time) |> 
      mutate(ecdf_value = (row_number()) / n()) |> 
      ungroup()
    
    #Create plot
    plot_speed_mode <- ggplot(ecdf_data, aes(x = last_time, y = ecdf_value, color = algorithm)) +
      geom_step(linewidth=1.6) +
      scale_color_manual(values = alg_colors, name = "Algorithm") +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(x = "seconds", y = "% of simulations")
    
    
    return(plot_speed_mode);
  }
  
  
}

##### Other files and functions to modify the data ##### 
{
  alg_correction <- tibble(alg=c("PT A-IIT(sq)","PT A-IIT(min)","PT_IIT_Z(sq)","PT_IIT_Z(min)","PT-IIT(sq)","PT-IIT(min)"),
                           algorithm=c("A-IIT","MH-mult","IIT","RF-MH","IIT","RF-MH"))
  fix_alg_name <- function(dataframe_input){
    dataframe_input <- dataframe_input |> 
      mutate(alg=str_remove(alg, "\\(\\d+\\)$")) |> #Modify name of algorithm, delete ID
      left_join(alg_correction,by="alg") #Add the standardized name of the algorithm
    return(dataframe_input)
  }
  
}

##### Create list for chosen_ids for each dimension #####
{
### Lowdim: First version of results
  chosen_ids_bimodal <- c(663,600,656,655)
  chosen_ids_7_mode <- c(614,612,615,613)
  
  chosen_dim <- "lowdim"
  lll <- create_plot_input("bimodal",chosen_ids_bimodal,chosen_dim)
  lll <- create_plot_input("7_mode",chosen_ids_7_mode,chosen_dim)
  
  mat_ids_ld <- rbind(chosen_ids_bimodal,chosen_ids_7_mode)
  list_ld_names <- c("bimodal","7_mode")
  
### Lowdim: Second version of lowdim
  chosen_ids_bimodal_new <- c(600,702,663,703)
  chosen_ids_7_mode_new <- c(612,700,614,701)
  
  chosen_dim <- "lowdim"
  lll <- create_plot_input("bimodal_new",chosen_ids_bimodal_new,chosen_dim)
  lll <- create_plot_input("7_mode_new",chosen_ids_7_mode_new,chosen_dim)
  
  
### Highdim: First version of results, only 25 simualtions
  chosen_ids_1k <- c(1040,1070,1042,1071)
  chosen_ids_3k <- c(936,937,938,939)
  chosen_ids_5k <- c(940,941,942,943)
  chosen_ids_7k <- c(944,945,946,947)
  
  chosen_dim <- "highdim"
  lll <- create_plot_input("dim_1k",chosen_ids_1k,chosen_dim)
  lll <- create_plot_input("dim_3k",chosen_ids_3k,chosen_dim)
  lll <- create_plot_input("dim_5k",chosen_ids_5k,chosen_dim)
  lll <- create_plot_input("dim_7k",chosen_ids_7k,chosen_dim)
  
  mat_ids_hd <- rbind(chosen_ids_1k,chosen_ids_3k,chosen_ids_5k,chosen_ids_7k)
  list_hd_names <- c(paste0("dim_",seq(1,7,by=2),"k"))
  
### Highdim: Second version of results, with 100 simulations
  chosen_ids_1k_100 <- c(1100:1103)
  chosen_ids_3k_100 <- c(1104:1107)
  chosen_ids_5k_100 <- c(1108:1111)
  chosen_ids_7k_100 <- c(1112:1115)
  
  chosen_dim <- "highdim"
  lll <- create_plot_input("dim_1k_100",chosen_ids_1k_100,chosen_dim)
  lll <- create_plot_input("dim_3k_100",chosen_ids_3k_100,chosen_dim)
  lll <- create_plot_input("dim_5k_100",chosen_ids_5k_100,chosen_dim)
  lll <- create_plot_input("dim_7k_100",chosen_ids_7k_100,chosen_dim)
  
}

##### Define the ids to use for the plots below
chosen_ids_1k <- c(1100:1103)
chosen_ids_3k <- c(1104:1107)
chosen_ids_5k <- c(1108:1111)
chosen_ids_7k <- c(1112:1115)

mat_ids_hd <- rbind(chosen_ids_1k,chosen_ids_3k,chosen_ids_5k,chosen_ids_7k)
mat_ids_hd_with_mult <-mat_ids_hd[,c(1,3)] #Ids ot algorithms that use multiplicity list
#used in Plots to compare algorithms that use Multiplicity list
list_hd_names <- c(paste0("dim_",seq(1,7,by=2),"k_100"))
dim_size <- c(1,3,5,7)

chosen_ids_bimodal <- c(600,702,663,703)
chosen_ids_7_mode <- c(612,700,614,701)

mat_ids_ld <- rbind(chosen_ids_bimodal,chosen_ids_7_mode)
list_ld_names <- c("bimodal_new","7_mode_new")
model_name <- c("bimodal","7_mode")

##### PLOTS and TABLES #####

#####  Plots to compare speed to modes of the 4 algorithms  #####
{
## Highdim
  chosen_dim <- "highdim"
  #Dim 1k
  j <- 1
  lll <- create_plot_input(list_hd_names[j],mat_ids_hd[j,],chosen_dim)
  (s_plot <- speed_plot_last(lll))
  export_plot(s_plot,"speed_to_mode",chosen_dim,dim_size[j])
  
  for(i in 1:4){
    lll <- create_plot_input(list_hd_names[i],mat_ids_hd[i,],chosen_dim)
    (s_plot <- speed_plot_last(lll))
    export_plot(s_plot,"speed_to_mode",chosen_dim,dim_size[i])
  }
  
#####################################################################  
## Lowdim
  chosen_dim <- "lowdim"
  lll <- create_plot_input(list_ld_names[1],mat_ids_ld[1,],chosen_dim)
  (s_plot <- speed_mode_lowdim(lll))
  export_plot(s_plot,paste0("speed_to_mode_",model_name[1]),chosen_dim,16)
  
  
  lll <- create_plot_input(list_ld_names[2],mat_ids_ld[2,],chosen_dim)
  (s_plot <- speed_mode_lowdim(lll))
  export_plot(s_plot,paste0("speed_to_mode_",model_name[2]),chosen_dim,16)

  # Other for bimodal
  lll <- create_plot_input("bimodal_other",c(600,663,704,705),chosen_dim)
  (s_plot <- speed_mode_lowdim(lll))
  
}

#### Plot of TVD ####
{
  chosen_dim <- "lowdim"
  #Low dim bimodal
  lll <- create_plot_input(list_ld_names[1],mat_ids_ld[1,],chosen_dim)
  (s_plot <- plot_tvd(lll,filter_measurement=50,filter_time=1000))
  export_plot(s_plot,"tvd_bimodal",chosen_dim,16)
  
  #Low dim 7_modes
  lll <- create_plot_input(list_ld_names[2],mat_ids_ld[2,],chosen_dim)
  (s_plot <- plot_tvd(lll,filter_measurement=50,filter_time=1000))
  export_plot(s_plot,"tvd_7_mode",chosen_dim,16)
  
  # Other for bimodal
  lll <- create_plot_input("bimodal_other",c(600,663,704,705),chosen_dim)
  (s_plot <- plot_tvd(lll,filter_measurement=50,filter_time=1000))
  
}


#####  Table: summary of swap_rate  #####
{
## For high dimension
  chosen_dim <- "highdim"
  
  for(i in 1:4){
    lll <- create_plot_input(list_hd_names[i],mat_ids_hd[i,])
    create_swap_rate_table(lll,chosen_dim,dim_size[i])
  }
  
###################################################################  
## For low dimension
 chosen_dim <- "lowdim"
 
 #Bimodal
 lll <- create_plot_input(list_ld_names[1],mat_ids_ld[1,],chosen_dim)
 create_swap_rate_table(lll,chosen_dim,dimension=model_name[1])
 # 7 modes
 lll <- create_plot_input(list_ld_names[2],mat_ids_ld[2,],chosen_dim)
 create_swap_rate_table(lll,chosen_dim,dimension=model_name[2])
 
 # Other for bimodal
 lll <- create_plot_input("bimodal_other",c(600,663,704,705),chosen_dim)
 # create_swap_rate_table(lll,chosen_dim,dimension="bimodal")
}

#####  Table: summary of iterations #####
{
# High dimension
  chosen_dim <- "highdim"
  
  for(i in 1:4){
    lll <- create_plot_input(list_hd_names[i],mat_ids_hd[i,],chosen_dim)
    create_table_iterations(lll,chosen_dim,dimension=dim_size[i])
  }

  
  ###################################################################
# Low dimension  
  chosen_dim <- "lowdim"
  for(i in 1:2){
    lll <- create_plot_input(list_ld_names[i],mat_ids_ld[i,],chosen_dim)
    create_table_iterations(lll,chosen_dim,dimension=model_name[i])
  }
  
  
}

#####  Table: report of temperatures #####
{
  
  # High dimension
  chosen_dim <- "highdim"
  for(i in 1:4){
    lll <- create_plot_input(list_hd_names[i],mat_ids_hd[i,],chosen_dim)
    create_table_temperatures(lll,chosen_dim,dimension=dim_size[i])
  }
  
  
  ###################################################################    
  # Low dimension  
  chosen_dim <- "lowdim"
  
  for(i in 1:2){
    lll <- create_plot_input(list_ld_names[i],mat_ids_ld[i,],chosen_dim)
    create_table_temperatures(lll,chosen_dim,dimension=model_name[i])
  }
  
  
}

#####  Plots to compare bounding constants  #####
{

  
  alg_correction <- tibble(alg=c("PT A-IIT(sq)","PT A-IIT(min)","PT_IIT_Z(sq)","PT_IIT_Z(min)"),
                           algorithm=c("A-IIT","MH-mult","IIT","RF-MH"))
  #Dim 3k
  chosen_ids <- c(932:935)+4
  lll <- create_plot_input("dim_3k",chosen_ids)
  
  bounds <- lll[["log_bounds"]]
  
  bounds <- bounds |> 
    mutate(bound=exp(log_bound),#Compute bounds from log_bounds
           alg=str_remove(alg, "\\(\\d+\\)$")) |> #Modify name of algorithm, delete ID
    left_join(alg_correction,by="alg") |> #Add the standardized name of the algorithm
    filter(algorithm=="A-IIT")# Consider only A-IIT beacuse MH doesnt have bounds
    
  #Beta labels for the x axis
  beta_labels <- sapply(sort(unique(bounds$replica)), function(i) {bquote(beta[.(i)])})
  
  bounds_plot <- bounds |> ggplot(aes(x=factor(replica),fill=algorithm, y=bound))+
    geom_violin()+
    scale_x_discrete(labels = beta_labels) +
    scale_fill_manual(values = alg_colors, name = "Algorithm") +
    labs(x = "", y = "bounding constant")+
    theme(legend.position = "none")
  
  
  export_plot(bounds_plot,"bounds_a_iit",3)
}


#####  Plots to compare algorithms that use Multiplicity list #####
{
  
  chosen_dim <- "highdim"
  for(i in 1:4){
    lll <- create_plot_input(list_hd_names[i],mat_ids_hd[i,])
    s_plot <- speed_plot_last(lll,mat_ids_hd_with_mult[i,])
    export_plot(s_plot,"speed_to_mode_RF",chosen_dim,dim_size[i])
  }
  

}

# These plots below consider previous IDs but they serve a different purpose, 
# just to compare the performance depending on number of RF replicas
# and the average number of iterations per replica
# Slight modifications in the temperatures don't affect that much these parameters

#####  Plots to compare number of RF replicas  #####
{
  #File with info to join
  add_info <- tibble(id=948:1007,
                     rf=rep(c(1:5,9:13),each=6),
                     dim=rep(c(3,3,5,5,7,7),10))
  add_info <- rbind(add_info,
                    tibble(id=1012:1023,
                           rf=rep(c(6,7),each=6),
                           dim=rep(c(3,3,5,5,7,7),2))
  )
  
  
### Create .Rds for the plots
  #dim 3k, 1-7,9-13
  chosen_ids <- c(seq(948,1007,by=6),
                  seq(949,1007,by=6),
                  seq(1012,1023,by=6),
                  seq(1013,1024,by=6))
  lll <- create_plot_input("dim_3k_rf_alg",chosen_ids)
  
  #dim 3k, 1-3
  # s_plot <- violin_plot_last(lll,add_info,subset_rfs=c(1:5))
  # s_plot
  # export_plot(s_plot,"compare_RF_1_5",3)

  
  #dim 3k compare RF side by side
  # ddd <- violin_plot_last(lll,add_info,return_dataset = T,subset_ids=chosen_ids)
  # 
  # s_plot <- ggplot(ddd,aes(x=factor(rf),fill=algorithm, y=time_visit))+
  #   geom_violin()+
  #   scale_fill_manual(values = alg_colors, name = "Algorithm")+
  #   labs(x = "Rejection Free replicas", y = "seconds")+
  #   facet_wrap(~algorithm)+
  #   theme(legend.position = "none",
  #         strip.text = element_text(size=25))
  # 
  # export_plot(s_plot,"compare_RF_side",dimension=3)
  
  #dim 3k compare RF for A-IIT
  ddd <- violin_plot_last(lll,add_info,return_dataset = T,subset_ids=chosen_ids)
  
  s_plot <- ddd |> filter(algorithm=="A-IIT") |> 
    ggplot(aes(x=factor(rf),fill=algorithm, y=time_visit))+
    geom_violin()+
    scale_fill_manual(values = alg_colors, name = "Algorithm")+
    labs(x = "Number of replicas using A-IIT", y = "seconds")+
    facet_wrap(~algorithm)+
    theme(legend.position = "none",
          strip.text = element_text(size=25))
  
  export_plot(s_plot,"compare_RF_a_iit","highdim",dimension=3)
  
  #dim 3k compare RF for A-IIT
  ddd <- violin_plot_last(lll,add_info,return_dataset = T,subset_ids=chosen_ids)
  
  s_plot <- ddd |> filter(algorithm=="MH-mult") |> 
    ggplot(aes(x=factor(rf),fill=algorithm, y=time_visit))+
    geom_violin()+
    scale_fill_manual(values = alg_colors, name = "Algorithm")+
    labs(x = "Number of replicas using MH with a multiplicity list", y = "seconds")+
    facet_wrap(~algorithm)+
    theme(legend.position = "none",
          strip.text = element_text(size=25))
  
  export_plot(s_plot,"compare_RF_mh_mult","highdim",dimension=3)
  
}

#####  Plot to compare number of iterations  #####
{
  chosen_id <- 1004 #dimension 5k
  lll <- create_plot_input("dim_5k_iterations",chosen_id)
  itera_data <- lll[["iterations"]]
  itera_data$replica <- as.numeric(itera_data$replica)
  algorithm_check <- str_extract(unique(itera_data$alg), "^[^(]+")
  if(algorithm_check=="PT A-IIT"){itera_data$alg <- "A-IIT"}else{cat("ERROR: not the correct algorithm")}
  
  itera_data |> ggplot(aes(x=as.factor(replica),y=iterations))+
    geom_violin()+
    scale_fill_manual(values = alg_colors, name = "Algorithm") +
    labs(x = "Rejection Free replicas", y = "seconds")
  
  
  beta_labels <- sapply(sort(unique(itera_data$replica)), function(i) {
    bquote(beta[.(i)])
  })
  
  itera_plot <- itera_data |> 
    ggplot(aes(x = as.factor(replica), y = iterations, fill=alg)) +
    geom_violin() +
    scale_x_discrete(labels = beta_labels) +
    scale_fill_manual(values = alg_colors, name = "Algorithm") +
    labs(x = "", y = "RF iterations")+
    theme(legend.position = "none",
          strip.text = element_text(size=25))
  
  export_plot(itera_plot,"compare_iterations")
  
  itera_data |> group_by(as.factor(replica)) |> summarise(avg.iter=mean(iterations))
  
}


#####  TITLE  #####
{
  
}
