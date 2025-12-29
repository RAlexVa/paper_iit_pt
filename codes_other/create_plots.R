rm(list=ls())
library(dplyr)
library(ggplot2)
library(stringr)
#### Create plots to include in the PDF
wd_path <- "C:/Users/ralex/Documents/src/paper_iit_pt"
source_path <- "C:/Users/ralex/Documents/src/paper_iit_pt/results"
export_path <- "C:/Users/ralex/Documents/src/paper-adaptive-iit-latex/images/highdim_ex"


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
create_list <- function(chosen_ids){
  ##### Define IDs and dimension #####
  chosen_dim <- "highdim"
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
    temperatures_temporal <- tibble(id=selected_id,temp_id=1:length(temperatures),temperature=temperatures);
    temperatures_matrix <- rbind(temperatures_matrix,temperatures_temporal)
    
    num_modes <- data_sum |> slice(i) |> pull(num_modes)
    #Add the BF to the algorithm name
    if(algorithm=="PT_A_IIT"){algorithm <- paste0("PT A-IIT(",data_sum |> slice(i)|> pull(bf),")")}
    if(algorithm=="PT_IIT_Z"){algorithm <- paste0("PT-IIT(",data_sum |> slice(i)|> pull(bf),")")}
    
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
    consolidated_list <- list()
    consolidated_list[["ids"]] <- chosen_ids;
    consolidated_list[["distances"]] <- distances;
    consolidated_list[["round_trip"]] <- round_trip;
    consolidated_list[["swap_rate"]] <- swap_rate;
    consolidated_list[["iterations"]] <- iterations;
    consolidated_list[["log_bounds"]] <- log_bounds;
    consolidated_list[["rf_replicas"]] <- rf_replicas;
    consolidated_list[["temperatures_matrix"]] <- temperatures_matrix;
    consolidated_list[["time_taken"]] <- time_taken;
    
    return(consolidated_list)

}
#Use the two functions above to import the list with information.
create_plot_input <- function(filename,input_ids){
  file_to_create <- file.path(source_path,paste0(filename,".Rds"))
  if(check_list(file_to_create,input_ids)){#If the list already exists and matches
    lista_creada <- readRDS(file_to_create);
  }else{#If the list doesn't exist or doesn't match
    lista_creada <- create_list(input_ids)
    saveRDS(create_list(input_ids),file_to_create);
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
    labs(x = "seconds", y = "%")
  
  
  
  
  
  
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

export_plot <- function(the_plot,plot_name,dimension){
  jpeg(file.path(export_path,paste0(plot_name,"_hd_",dimension,"k.jpg")),
       width=1200,height = 600,pointsize = 30)
  print(the_plot)
  dev.off()
  message("Successfully created plot: ",paste0(plot_name,"_hd_",dimension,"k.jpg"))
}


##### PLOTS #####

#####  Plots to compare 4 algorithms  #####
{
  #Dim 1k
  chosen_ids <- c(932:935)
  lll <- create_plot_input("dim_1k",chosen_ids)
  s_plot <- speed_plot_last(lll)
  export_plot(s_plot,"speed_to_mode",1)
  #Dim 3k
  chosen_ids <- c(932:935)+4
  lll <- create_plot_input("dim_3k",chosen_ids)
  s_plot <- speed_plot_last(lll)
  export_plot(s_plot,"speed_to_mode",3)
  #Dim 5k
  chosen_ids <- c(932:935)+8
  lll <- create_plot_input("dim_5k",chosen_ids)
  s_plot <- speed_plot_last(lll)
  export_plot(s_plot,"speed_to_mode",5)
  #Dim 7k
  chosen_ids <- c(932:935)+12
  lll <- create_plot_input("dim_7k",chosen_ids)
  s_plot <- speed_plot_last(lll)
  export_plot(s_plot,"speed_to_mode",7)
  
}


#####  Plots to compare RF algorithms #####
{
  #Dim 1k
  chosen_ids <- c(932:935)
  lll <- create_plot_input("dim_1k",chosen_ids)
  s_plot <- speed_plot_last(lll,c(932,934))
  export_plot(s_plot,"speed_to_mode_RF",1)
  #Dim 3k
  chosen_ids <- c(932:935)+4
  lll <- create_plot_input("dim_3k",chosen_ids)
  s_plot <- speed_plot_last(lll,c(932,934)+4)
  export_plot(s_plot,"speed_to_mode_RF",3)
  #Dim 5k
  chosen_ids <- c(932:935)+8
  lll <- create_plot_input("dim_5k",chosen_ids)
  s_plot <- speed_plot_last(lll,c(932,934)+8)
  export_plot(s_plot,"speed_to_mode_RF",5)
  #Dim 7k
  chosen_ids <- c(932:935)+12
  lll <- create_plot_input("dim_7k",chosen_ids)
  s_plot <- speed_plot_last(lll,c(932,934)+12)
  export_plot(s_plot,"speed_to_mode_RF",7)
}


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
  s_plot <- violin_plot_last(lll,add_info,subset_rfs=c(1:5))
  s_plot
  export_plot(s_plot,"compare_RF_1_5",3)

  
  #dim 3k compare RF side by side
  ddd <- violin_plot_last(lll,add_info,return_dataset = T,subset_ids=chosen_ids)
  
  s_plot <- ggplot(ddd,aes(x=factor(rf),fill=algorithm, y=time_visit))+
    geom_violin()+
    scale_fill_manual(values = alg_colors, name = "Algorithm")+
    labs(x = "Rejection Free replicas", y = "seconds")+
    facet_wrap(~algorithm)+
    theme(legend.position = "none",
          strip.text = element_text(size=25))
  
  export_plot(s_plot,"compare_RF_side",dimension=3)
  
  #dim 3k compare RF for A-IIT
  ddd <- violin_plot_last(lll,add_info,return_dataset = T,subset_ids=chosen_ids)
  
  s_plot <- ddd |> filter(algorithm=="A-IIT") |> 
    ggplot(aes(x=factor(rf),fill=algorithm, y=time_visit))+
    geom_violin()+
    scale_fill_manual(values = alg_colors, name = "Algorithm")+
    labs(x = "Rejection Free replicas", y = "seconds")+
    facet_wrap(~algorithm)+
    theme(legend.position = "none",
          strip.text = element_text(size=25))
  
  export_plot(s_plot,"compare_RF_a_iit",dimension=3)
  
  
  
}

#####  TITLE  #####
{
  
}

#####  TITLE  #####
{
  
}

#####  TITLE  #####
{
  
}
