##### Start of the code #####
# This code works with output from the highd_seeded.R
#Each .Rds file contains results for 1 simulation
#In the name of the file it's the ID from the input file and the seed for the simulation
# Do ctrl+F IMPORTANT to find how using this code with a different output might break the code

rm(list=ls())
# Import libraries
library(stringr)
library(tidyverse)
library(gridExtra)# For tables
library(latex2exp) #For using latex
library(survminer); library(survival) #For plot showing how fast each replica visits all modes



##### Define IDs and dimension #####
# Choose dimension
# chosen_dim <- "highdim"; file_dim <- "highd"
chosen_dim <- "lowdim"; file_dim <- "lowd"
chosen_ids <- c(600:615,616:617,620:635,655:662,663:668)
# chosen_ids <- c(1070:1071)#c(932:1007)
# chosen_ids <- c(908:911)#Spaced model
##### Read files and specifications#####
#Read CSV with simulation details
parameters <- read_csv(paste0("inputs/simulation_details_",file_dim,".csv"), col_types = cols())

#List of files in the results folder
data_sum <- tibble(file_names=list.files(path = "results", pattern = "^sim_.*\\.Rds")) |> 
  mutate(id=as.numeric(str_extract(file_names, "id_([0-9]+)",group=1)),
         sim_id=str_extract(file_names, "id_[0-9]+_([0-9]+)\\.Rds",group=1),
         dim=str_extract(file_names, "(?<=sim_)[^_]+(?=_id)")) |> 
  filter(dim==chosen_dim) |> 
  left_join(parameters, by="id")|> 
  filter(id %in% chosen_ids)


##### Create dataset to store results#####
if(chosen_dim=="highdim"){
  # Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp") #To use eval_lik function
  # source("functions/r_functions.R") #To get dimension of the problem
  
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
}
if(chosen_dim=="lowdim"){
  
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
  # rf_replicas <- tibble(id=numeric(),rf_reps=numeric());
  temperatures_matrix <- tibble(id=numeric(),temp_id=numeric(),temperature=numeric());
  
  #These are only for lowdim results
  time_visit <- tibble(alg=character(),sim=numeric(),mode=character(),time=numeric())
  #time visit measures time for beta=1 to visit each state
  tvd_report <- tibble(alg=character(),sim=numeric(),measurement=numeric(),tvd=numeric(),time=numeric());
  #Extimated probability for the modes
  prob_modes <- tibble(alg=character(),sim=numeric(),mode=character(),prob=numeric());
  
  rf_replicas <- tibble(id=numeric(),rf_reps=numeric());
  #Round trips: Number of roundtrips of each replica
  #tvd: Final TVD  
  #pi_modes: estimated density of the modes
  #mode_visit: Iteration of the first visit to the mode

  mod1 = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1);
  mod2 = c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0);
  mod3 = c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1);
  mod4 = c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0);
  mod5 = c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1);
  mod6 = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1);
  mod7 = c(0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0);
  
  matrix_mode <- cbind(mod1,mod2,mod3,mod4,mod5,mod6,mod7)
  index_mode <- c(65536,21846,43691,256,65281,32770,385)
  
}
time_taken <- as.data.frame(matrix(nrow=0,ncol=2));colnames(time_taken) <- c("alg","time")

##### FOR loop to read .Rds files and get results#####
# Start creating datasets with information
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
    # number_replicas_print <- temp_total_number
    # temperatures_to_print <- paste(temp_total_number,
    #                                " temperatures. From: ",
    #                                temp_from,
    #                                " To: ",
    #                                temp_to,collapse=',')
  }
  
  temperatures_temporal <- tibble(id=selected_id,temp_id=1:length(temperatures),temperature=temperatures);
  temperatures_matrix <- rbind(temperatures_matrix,temperatures_temporal)
  
  if(chosen_dim=="lowdim"){
    check_model <- data_sum |> slice(i) |> pull(model)
    if(check_model=="bimodal"){num_modes <- 2; spec_matrix_mode <- matrix_mode[,2:3];spec_index_mode <- index_mode[2:3];}
    if(check_model=="7_mode"){num_modes <- 7;spec_matrix_mode <- matrix_mode;spec_index_mode <- index_mode;}
  }
  if(chosen_dim=="highdim"){
    num_modes <- data_sum |> slice(i) |> pull(num_modes)
  }
  
  
  #Add the BF to the algorithm name
  if(algorithm=="PT_A_IIT"){algorithm <- paste0("PT A-IIT(",data_sum |> slice(i)|> pull(bf),")")}
  if(algorithm=="PT_IIT_Z"){algorithm <- paste0("PT-IIT(",data_sum |> slice(i)|> pull(bf),")")}
  
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

}


##### Export plots and tables #####
export_path <- paste0("C:/Users/ralex/Documents/src/paper_iit_pt/images/",chosen_dim,"_ex")
export_file_name <- paste0(paste0(chosen_ids,collapse="_"),"_",chosen_dim)

#Starts high dim reports
if(chosen_dim=="highdim"){
  
  
##### Add some variables (id, algorithm)
  distances$id <- str_extract(distances$alg, "\\(\\d+\\)") |> 
    str_remove_all("[()]")  # Extract numbers between parentheses and remove the parentheses
  distances$algorithm <- str_replace(distances$alg, "\\s*\\(.*", "")  # Extract text before parenthesis
  
##### REPORT: Speed to modes for maximum temperature #####
  {
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
    # mutate(first_time=min(m1,m2),last_time=max(m1,m2))
    dist_t1_times[dist_t1_times==Inf] <- NA ## replace Infinities for NAs
    #How many didn't reach the optimum
    table(dist_t1$min_dist)

    
    
##### REPORT: First for the highest inverse temperature to visit all modes #####
    {
      report_first <- dist_t1_times |> 
        select(alg,first_time) |> 
        rename(time_visit=first_time)
      
      report_first_bk <- report_first #Backup
      ids_to_print <- chosen_ids
      ids_to_print <- c(932:935)+12#+4#+8
      ids_to_print <- c(932,934)+12#+4#+8
      ids_to_print <- seq(948,1007,by=6)
      report_first <- report_first_bk |> filter(grepl(paste(ids_to_print,collapse = "|"),alg))
      fit <- survfit(Surv(time_visit,rep(1,nrow(report_first)))~alg,data=report_first)
      
      
      (plot_surv_mode <- ggsurvplot(fit,
                                    data=report_first,
                                    fun="event",
                                    palette = "Set1",    # Color palette
                                    xlab = "Time (seconds)",
                                    ylab = "Prop. of simulations visiting all modes",
                                    legend.title = "Algorithm",
                                    # break.time.by = time_br,
                                    font.x = 25,        # X-axis label font size
                                    font.y = 25,        # Y-axis label font size
                                    font.tickslab = 18, # Axis tick labels (numbers) font size
                                    font.legend = 20,
                                    size=2,
                                    conf.int = FALSE,
                                    censor = TRUE))   # Legend text font size)
      
      export_file_name_temp <- export_file_name
      if(!identical(ids_to_print,chosen_ids)){
        export_file_name_temp <- paste0(paste0(ids_to_print,collapse="_"),"_",chosen_dim)}
      jpeg(file.path(export_path,paste0(export_file_name_temp,"_speed_mode_1st_rep_1st_visit",".jpg")),width=1200,height =600,pointsize = 30)
      print(plot_surv_mode)
      dev.off()
      
}
    
##### REPORT: Last time for the highest inverse temperature to visit all modes #####    
    {
      report_last <- dist_t1_times |> 
        select(alg,last_time)|> 
        rename(time_visit=last_time)
      
      report_last_bk <- report_last
      
      
      ids_to_print <- chosen_ids
      ids_to_print <- c(1024:1027)+8#+12#+4#+8
      ids_to_print <- c(1024,1026)#+12#+4#+8
      report_last <- report_last_bk |> filter(grepl(paste(ids_to_print,collapse = "|"),alg)) |> 
        filter(time_visit<100000)
      fit <- survfit(Surv(time_visit,rep(1,nrow(report_last)))~alg,data=report_last)
      
      
      (plot_surv_mode <- ggsurvplot(fit,
                                    data=report_last,
                                    fun="event",
                                    palette = "Set1",    # Color palette
                                    xlab = "Time (seconds)",
                                    ylab = "Prop. of simulations visiting all modes",
                                    legend.title = "Algorithm",
                                    # break.time.by = time_br,
                                    font.x = 25,        # X-axis label font size
                                    font.y = 25,        # Y-axis label font size
                                    font.tickslab = 18, # Axis tick labels (numbers) font size
                                    font.legend = 20,
                                    size=2,
                                    conf.int = FALSE,
                                    censor = TRUE))   # Legend text font size)
      
      export_file_name_temp <- export_file_name
      if(!identical(ids_to_print,chosen_ids)){
        export_file_name_temp <- paste0(paste0(ids_to_print,collapse="_"),"_",chosen_dim)}
      jpeg(file.path(export_path,paste0(export_file_name_temp,"_speed_mode_1st_rep_last_visit",".jpg")),width=1200,height =600,pointsize = 30)
      print(plot_surv_mode)
      dev.off()
      
    }
    
    
##### REPORT: First for the highest inverse temperature COMPARING RF replicas #####
    {
      report_first_compare <- dist_t1_times |> 
        select(alg,id,first_time) |> 
        rename(time_visit=first_time) |> 
        mutate(id=as.numeric(id)) 
      
      report_first_compare$bf <- str_extract(report_first_compare$alg, "(?<=\\()[^)]+(?=\\))")
      
      #Filter the IDs to compare by number of RF replicas
      report_first_compare <- report_first_compare |> filter(id>=948)
      
#Data set with additional information
      add_info <- tibble(id=948:1007,
                         rf=rep(c(1:5,9:13),each=6),
                         dim=rep(c(3,3,5,5,7,7),10))
      
      
 #Dataset with all the info     
      report_firstc_bk <- report_first_compare |>  #Backup
                          left_join(add_info,by="id")
      
# Parameters to filter      
      chosen_bf <- "min"; #chosen_bf <- "min";
      dim_size <- 7#3#5#7
      report_first_compare <- report_firstc_bk |> 
        filter(bf==chosen_bf) |> 
        filter(dim==dim_size)
      
      fit <- survfit(Surv(time_visit,rep(1,nrow(report_first_compare)))~id+rf,data=report_first_compare)
      
      
      (plot_surv_mode <- ggsurvplot(fit,
                                    data=report_first_compare,
                                    fun="event",
                                    palette = "Set3",    # Color palette
                                    xlab = "Time (seconds)",
                                    ylab = "Prop. of simulations visiting all modes",
                                    legend.title = "Alg",
                                    # break.time.by = time_br,
                                    font.x = 25,        # X-axis label font size
                                    font.y = 25,        # Y-axis label font size
                                    font.tickslab = 18, # Axis tick labels (numbers) font size
                                    font.legend = 20,
                                    size=2,
                                    conf.int = FALSE,
                                    censor = TRUE))   # Legend text font size)
      
      
      export_file_name_temp <- paste0("first_visit_compare_rf_",chosen_bf,"_",dim_size,".jpg")
      
      jpeg(file.path(export_path,export_file_name_temp),width=1200,height =600,pointsize = 30)
      print(plot_surv_mode)
      dev.off()
      
    }    

##### REPORT: Last for the highest inverse temperature COMPARING RF replicas #####
    {
      report_first_compare <- dist_t1_times |> 
        select(alg,id,last_time) |> 
        rename(time_visit=last_time) |> 
        mutate(id=as.numeric(id)) 
      
      report_first_compare$bf <- str_extract(report_first_compare$alg, "(?<=\\()[^)]+(?=\\))")
      
      #Filter the IDs to compare by number of RF replicas
      report_first_compare <- report_first_compare |> filter(id>=948)
      
      #Data set with additional information
      add_info <- tibble(id=948:1007,
                         rf=rep(c(1:5,9:13),each=6),
                         dim=rep(c(3,3,5,5,7,7),10))
      
      
      #Dataset with all the info     
      report_firstc_bk <- report_first_compare |>  #Backup
        left_join(add_info,by="id")
      
      # Parameters to filter      
      chosen_bf <- "min"; #chosen_bf <- "min";
      dim_size <- 7#3#5#7
      report_first_compare <- report_firstc_bk |> 
        filter(bf==chosen_bf) |> 
        filter(dim==dim_size)
      
      fit <- survfit(Surv(time_visit,rep(1,nrow(report_first_compare)))~id+rf,data=report_first_compare)
      
      
      (plot_surv_mode <- ggsurvplot(fit,
                                    data=report_first_compare,
                                    fun="event",
                                    palette = "Set3",    # Color palette
                                    xlab = "Time (seconds)",
                                    ylab = "Prop. of simulations visiting all modes",
                                    legend.title = "Alg",
                                    # break.time.by = time_br,
                                    font.x = 25,        # X-axis label font size
                                    font.y = 25,        # Y-axis label font size
                                    font.tickslab = 18, # Axis tick labels (numbers) font size
                                    font.legend = 20,
                                    size=2,
                                    conf.int = FALSE,
                                    censor = TRUE))   # Legend text font size)
      
      
      export_file_name_temp <- paste0("last_visit_compare_rf_",chosen_bf,"_",dim_size,".jpg")
      
      jpeg(file.path(export_path,export_file_name_temp),width=1200,height =600,pointsize = 30)
      print(plot_surv_mode)
      dev.off()
      
    }           

##### REPORT: Last for the highest inverse temperature COMPARING RF replicas AND BF #####
    {
      report_first_compare <- dist_t1_times |> 
        select(alg,id,last_time) |> 
        rename(time_visit=last_time) |> 
        mutate(id=as.numeric(id)) 
      
      report_first_compare$bf <- str_extract(report_first_compare$alg, "(?<=\\()[^)]+(?=\\))")
      
      #Filter the IDs to compare by number of RF replicas
      report_first_compare <- report_first_compare |> filter(id>=948)
      
      #Data set with additional information
      add_info <- tibble(id=948:1007,
                         rf=rep(c(1:5,9:13),each=6),
                         dim=rep(c(3,3,5,5,7,7),10))
      
      
      #Dataset with all the info     
      report_firstc_bk <- report_first_compare |>  #Backup
        left_join(add_info,by="id")
      
      # Parameters to filter      
      dim_size <- 7#3#5#7
      rf_select <- c(10:13)
      report_first_compare <- report_firstc_bk |> 
        filter(dim==dim_size) |> 
        filter(rf%in%rf_select)
      
      fit <- survfit(Surv(time_visit,rep(1,nrow(report_first_compare)))~bf+rf,data=report_first_compare)
      
      
      (plot_surv_mode <- ggsurvplot(fit,
                                    data=report_first_compare,
                                    fun="event",
                                    palette = "Set3",    # Color palette
                                    xlab = "Time (seconds)",
                                    ylab = "Prop. of simulations visiting all modes",
                                    legend.title = "Alg",
                                    # break.time.by = time_br,
                                    font.x = 25,        # X-axis label font size
                                    font.y = 25,        # Y-axis label font size
                                    font.tickslab = 18, # Axis tick labels (numbers) font size
                                    font.legend = 20,
                                    size=2,
                                    conf.int = FALSE,
                                    censor = TRUE))   # Legend text font size)
      
      
      export_file_name_temp <- paste0("last_visit_compare_BFS_",dim_size,"_",paste0(rf_select,collapse="_"),".jpg")
      
      jpeg(file.path(export_path,export_file_name_temp),width=1200,height =600,pointsize = 30)
      print(plot_surv_mode)
      dev.off()
      
    }         
    
        
  }
  
##### REPORT: minimum distance reached by ANY replica in ANY simulation #####
#Used in case there are some simulations that didn't finish.  
  {
    report_min_dist <- distances |> 
      group_by(alg,mode) |> 
      slice_min(min_dist,n=1) |> 
      slice_min(time_find,n=1) |> 
      select(alg,mode,min_dist,time_find) |> 
      ungroup()
    # view(report_min_dist)   
    
    wide_report_min_dist <-  report_min_dist |> 
      select(-time_find) |> 
      pivot_wider(names_from = mode,values_from = min_dist)
    view(wide_report_min_dist)
    
    
    
    report_by_alg_sim <- distances |> 
      group_by(alg,mode,sim) |> 
      slice_min(min_dist,n=1) |> 
      slice_min(time_find,n=1) |> 
      select(alg,mode,min_dist,time_find) |> 
      ungroup() |> 
      select(-time_find) |> 
      pivot_wider(names_from = mode,values_from = min_dist)
    
    
    modes_not_visited <- report_by_alg_sim |> 
      rowwise() |> 
      mutate(max_dist=max(across(starts_with("m")),na.rm=T)) |> 
      filter(max_dist>0)
    
    #Export table with averages
    jpeg(file.path(export_path,paste0("sim_not_find_modes",".jpg")),width=80*ncol(modes_not_visited),height =26*nrow(modes_not_visited),pointsize = 30)
    grid.arrange(tableGrob(modes_not_visited))
    dev.off()  
    
    
    
  }
   
  
##### REPORT: Time of first visit of modes by ANY replica #####
  #For each simulation the time of first visit by any replica and the replica that did it.
  {
    dist_mode_times <- distances  |> 
      mutate(time_fixed=ifelse(min_dist>0,Inf,time_find)) |> 
      select(-min_dist,-time_find) |> 
      pivot_longer(cols=c(time_fixed),names_to="variable",values_to="measure") |> 
      pivot_wider(names_from=mode,values_from=measure)|> 
      rowwise() |> 
      mutate(min_time = min(across(matches("^m\\d+$")))) |> 
      filter(min_time<Inf)
    
    first_visit_report <- dist_mode_times |> 
      select(-variable,-min_time) |> 
      group_by(alg,algorithm,id,sim) |> 
      summarise(across(starts_with("m"),~min(.x,na.rm=T),.names="min_{.col}"),
                across(starts_with("m"),~ {
                  # col_name <- cur_column()
                  min_val <- min(.x, na.rm = TRUE)
                  # Find the first row where this minimum value occurs
                  first(which(.x == min_val))
                },.names="row_{.col}")
      ) |> 
      select(-starts_with("row_min_")) |> 
      rowwise() |> 
      mutate(last_visit = max(across(starts_with("min_m")))) |> 
      ungroup()
    
    ##### REPORT: Number of modes visited in each simulation #####  
    {
      report_modes_visited <- first_visit_report |> 
        select(alg,algorithm,id,sim,starts_with("min_m")) |> 
        pivot_longer(cols=starts_with("min_m"),names_to = "mode",values_to="time") |> 
        mutate(visited=time!=Inf) |> 
        group_by(id,alg,sim) |> 
        summarise(modes_visited=sum(visited)) |> 
        ungroup() 
      
      (summary_modes_visited <-   report_modes_visited |> 
          group_by(alg,modes_visited) |> 
          summarise(count_sims=n()))
      ##### REPORT: Id's that didn't visit a specific number of modes #####  
      tot_modes <- 6
      sim_not_visiting <- report_modes_visited |> 
        filter(modes_visited!=tot_modes) |> 
        pull(sim)
      View(sim_not_visiting)
    }    
    
  
    
    ##### REPORT: Speed to mode considering that any replica visits the mode
    #Minimum time to visit all the modes, the modes can be visited by any replica
      #Report of speed to mode considering that any replica visits the mode
      forsurv <- first_visit_report |> select(alg,last_visit)
      forsurv <- forsurv |> filter(last_visit<Inf) #Filter only the ones that finished
     
      forsurv_bk <- forsurv #Backup
      ids_to_print <- chosen_ids
      ids_to_print <- c(884,886)+8#c(870,872)#c(880,881,806,808)
      forsurv <- forsurv_bk |> filter(grepl(paste(ids_to_print,collapse = "|"),alg))
      fit <- survfit(Surv(last_visit,rep(1,nrow(forsurv)))~alg,data=forsurv)
      
      (plot_surv_mode <- ggsurvplot(fit,
                                    data=forsurv,
                                    fun="event",
                                    palette = "Set1",    # Color palette
                                    xlab = "Time (seconds)",
                                    ylab = "Prop. of simulations visiting all modes",
                                    legend.title = "Algorithm",
                                    # break.time.by = time_br,
                                    font.x = 25,        # X-axis label font size
                                    font.y = 25,        # Y-axis label font size
                                    font.tickslab = 18, # Axis tick labels (numbers) font size
                                    font.legend = 20,
                                    size=2,
                                    conf.int = FALSE,
                                    censor = TRUE))   # Legend text font size)
      
      export_file_name_temp <- export_file_name
      if(!identical(ids_to_print,chosen_ids)){
        export_file_name_temp <- paste0(paste0(ids_to_print,collapse="_"),"_",chosen_dim)}
      jpeg(file.path(export_path,paste0(export_file_name_temp,"_speed_mode_anyrep",".jpg")),width=1200,height =600,pointsize = 30)
      print(plot_surv_mode)
      dev.off()
      }
    
  #Modify the name of the file to export
  if(nchar(export_file_name)>80){
    export_file_name <- paste0(min(chosen_ids),"_",max(chosen_ids),"_many_ids_",chosen_dim)
  }

   
  
  
  
  
}
  

if(chosen_dim=="lowdim"){
  ##### REPORT: Time to visit modes
  time_visit$id <- as.numeric(str_extract(time_visit$alg, "\\(\\d+\\)") |> 
                               str_remove_all("[()]"))
  ###Check if one mode was not visited
  time_visit |> filter(time==0);
  time_visit |> filter(is.na(time));
  # Now start the report
  time_visit_report <- time_visit |> 
    pivot_wider(names_from=mode, values_from=time) |> 
    mutate(count_modes=rowSums(!is.na(across(starts_with("m"))))) |> 
    rowwise() |> 
    mutate(first_time = min(across(matches("^m\\d+$")),na.rm=T),
           last_time = max(across(matches("^m\\d+$")),na.rm=T))
  
  # Just to check that all the modes (either 2 or 7) were visited
  time_visit_report |> group_by(count_modes) |> summarise(total_count=n())
  
 ##### REPORT: Only for bimodal 
  report_bimodal <- time_visit_report |> 
    filter(count_modes==2) |> 
    select(id,alg,last_time)

  report_bk <- report_bimodal
  ids_to_print <- c(663,665) #MH mult bimodal
  ids_to_print <- c(600) #A-IIT (sq) bimodal
  ids_to_print <- c(656) #RF-MH bimodal
  ids_to_print <- c(655) #PT-IIT bimodal
  ids_to_print <- c(663,600,656,655) #Final chosen for bimodal
  report_bk <- report_bimodal |> filter(id %in% ids_to_print)
    
  fit <- survfit(Surv(last_time,rep(1,nrow(report_bk)))~alg,data=report_bk)

  
  (plot_surv_mode <- ggsurvplot(fit,
                                data=report_bk,
                                fun="event",
                                palette = "Set1",    # Color palette
                                xlab = "Time (seconds)",
                                ylab = "Prop. of simulations visiting all modes",
                                legend.title = "Algorithm",
                                # break.time.by = time_br,
                                font.x = 25,        # X-axis label font size
                                font.y = 25,        # Y-axis label font size
                                font.tickslab = 18, # Axis tick labels (numbers) font size
                                font.legend = 20,
                                size=2,
                                conf.int = FALSE,
                                censor = TRUE))   # Legend text font size)
  
  export_file_name_temp <- export_file_name
  if(!identical(ids_to_print,chosen_ids)){
    export_file_name_temp <- paste0(paste0(ids_to_print,collapse="_"),"_",chosen_dim)}


  jpeg(file.path(export_path,paste0(export_file_name_temp,"_speed_mode_bimodal",".jpg")),width=1200,height =600,pointsize = 30)
  print(plot_surv_mode)
  dev.off()
  
  ##### REPORT: Only for 7-modes 
  report_multimodal <- time_visit_report |> 
    filter(count_modes==7) |> 
    select(id,alg,last_time)
  
  report_bk <- report_multimodal
  ids_to_print <- c(606,614) #MH-Mult 7-mode
  ids_to_print <- c(612) #PT A-IIT 7-mode
  ids_to_print <- c(615) #RF-MH 7-mode
  ids_to_print <- c(613) #PT-IIT 7-mode
  ids_to_print <- c(614,612,615,613) #Final chosen for 7-mode
  report_bk <- report_multimodal |> filter(id %in% ids_to_print)
  
  
  fit <- survfit(Surv(last_time,rep(1,nrow(report_bk)))~alg,data=report_bk)
  
  
  (plot_surv_mode <- ggsurvplot(fit,
                                data=report_bk,
                                fun="event",
                                palette = "Set1",    # Color palette
                                xlab = "Time (seconds)",
                                ylab = "Prop. of simulations visiting all modes",
                                legend.title = "Algorithm",
                                # break.time.by = time_br,
                                font.x = 25,        # X-axis label font size
                                font.y = 25,        # Y-axis label font size
                                font.tickslab = 18, # Axis tick labels (numbers) font size
                                font.legend = 20,
                                size=2,
                                conf.int = FALSE,
                                censor = TRUE))   # Legend text font size)
  
  export_file_name_temp <- export_file_name
  if(!identical(ids_to_print,chosen_ids)){
    export_file_name_temp <- paste0(paste0(ids_to_print,collapse="_"),"_",chosen_dim)}
  
  jpeg(file.path(export_path,paste0(export_file_name_temp,"_speed_mode_7mode",".jpg")),width=1200,height =600,pointsize = 30)
  print(plot_surv_mode)
  dev.off()
  
  
  
  
  
  

  
  #### REPORT: TVD vs time
  
  # For each measurement we average the time and average the TVD
  tvd_report$id <- as.numeric(str_extract(tvd_report$alg, "\\(\\d+\\)") |> 
                                str_remove_all("[()]"))
  tvd_summary <- tvd_report |> 
    group_by(id,alg,measurement) |> 
    summarise(avg_tvd=mean(tvd),avg_time=mean(time)) |> 
    ungroup()
  
  
  ids_to_print <- c(663,665) #MH mult bimodal
  ids_to_print <- c(600) #A-IIT (sq) bimodal
  ids_to_print <- c(656) #RF-MH bimodal
  ids_to_print <- c(655) #PT-IIT bimodal
  ids_to_print <- c(663,600,656,655) #Final chosen for bimodal
  
  ids_to_print <- c(606,614) #MH-Mult 7-mode
  ids_to_print <- c(612) #PT A-IIT 7-mode
  ids_to_print <- c(615) #RF-MH 7-mode
  ids_to_print <- c(613) #PT-IIT 7-mode
  ids_to_print <- c(614,612,615,613) #Final chosen for 7-mode
  
  (tvd_plot <- tvd_summary |> filter(id %in% ids_to_print) |> 
    filter(measurement>50) |>
    ggplot(aes(x=avg_time,y=avg_tvd,color=alg))+
    geom_line())
  
  export_file_name_temp <- export_file_name
  if(!identical(ids_to_print,chosen_ids)){
    export_file_name_temp <- paste0(paste0(ids_to_print,collapse="_"),"_",chosen_dim)}
  
  
  jpeg(file.path(export_path,paste0(export_file_name_temp,"_tvd_time",".jpg")),width=1200,height =600,pointsize = 30)
  print(tvd_plot)
  dev.off()
  

  
}

export_file_name <- "6XX_many_ids_lowdim"

##### REPORT: Average number of iterations per algorithm #####
{
  iterations$id <- as.numeric(str_extract(iterations$alg, "\\(\\d+\\)") |> 
                                str_remove_all("[()]"))
  iterations$replica <- as.numeric(iterations$replica)
  
  
  iter_report <- iterations |> 
    mutate(replica=as.numeric(replica)) |> 
    group_by(alg,id,replica) |> 
    summarise(avg.iter=mean(iterations,na.rm=T)) |> ungroup() |> 
    arrange(id) |> select(-id)
  
  #Horizontal
  horizontal_report_iter <- iter_report |> 
    pivot_wider(names_from = replica,values_from = avg.iter)
  # View(horizontal_report_iter)
  
  #Vertical
  vertical_report_iter <- iter_report |> 
    pivot_wider(names_from = alg,values_from = avg.iter)
  # View(vertical_report_iter)
  
  #Choose which report to export
  table_iter_report <- vertical_report_iter
  # table_iter_report <- horizontal_report_iter
  
  iterations$id <- as.numeric(str_extract(iterations$alg, "\\(\\d+\\)") |> 
                                str_remove_all("[()]"))  #Extract ID
  
  
  boxplot_iter <- iterations |> 
    left_join(rf_replicas, by="id") |> #Add the number of Rejection Free temperatures
    mutate(replica=as.numeric(replica)) |> #Transform replica id into number
    filter(replica<=rf_reps) |> #Filter only the iterations for RF replicas
    ggplot(aes(x=factor(replica),y=iterations))+
    geom_boxplot()+
    facet_wrap(~alg, scales="free")
  # facet_wrap(~alg)
  boxplot_iter
  
  #Export table with averages
  jpeg(file.path(export_path,paste0(export_file_name,"_iterations",".jpg")),width=140*ncol(table_iter_report),height =28*nrow(table_iter_report),pointsize = 30)
  grid.arrange(tableGrob(table_iter_report |> mutate(across(starts_with("P"),\(x) round(x,2)))))
  dev.off()  
  
  # Export boxplot
  jpeg(file.path(export_path,paste0(export_file_name,"_iter_boxplot",".jpg")),width=1200,height =600,pointsize = 30)
  print(boxplot_iter)
  dev.off()  
  
  
}   

##### REPORT: Average swap rate #####
{
  swap_rate$id <- as.numeric(str_extract(swap_rate$alg, "\\(\\d+\\)") |> 
                               str_remove_all("[()]"))
  swap_rate$replica <- as.numeric(swap_rate$replica)
  
  
  sr_report <- swap_rate |>
    group_by(id,alg,replica) |> 
    summarise(avg.sr=mean(swap_rate,na.rm=T)) |> ungroup() |> 
    arrange(id,replica)
  
  #Horizontal
  horizontal_report_sr <- sr_report |> select(-id) |> 
    pivot_wider(names_from = replica,values_from = avg.sr)
  # View(horizontal_report_sr)
  
  #Vertical
  vertical_report_sr <- sr_report |> select(-id) |> 
    pivot_wider(names_from = alg,values_from = avg.sr)
  # View(vertical_report_sr)
  
  
  boxplot_sr <- swap_rate |> 
    ggplot(aes(x=factor(replica),y=swap_rate))+
    geom_boxplot()+
    geom_hline(yintercept=0.234,col="red")+
    facet_wrap(~alg)
  boxplot_sr
  
  table_sr_report <- vertical_report_sr
  # table_sr_report <- horizontal_report_sr
  # View(table_sr_report)
  
  
  #Export table with averages
  jpeg(file.path(export_path,paste0(export_file_name,"_swap_rate",".jpg")),width=120*ncol(table_sr_report),height =30*nrow(table_sr_report),pointsize = 30)
  grid.arrange(tableGrob(table_sr_report |> mutate(across(starts_with("P"),\(x) round(x,6)))))
  dev.off()  
  
  # Export boxplot
  jpeg(file.path(export_path,paste0(export_file_name,"_sr_boxplot",".jpg")),width=1200,height =600,pointsize = 30)
  print(boxplot_sr)
  dev.off()  
  
  
  
}  


##### REPORT: Round trips#####
{
  
  round_trip$id <- as.numeric(str_extract(round_trip$alg, "\\(\\d+\\)") |> 
                                str_remove_all("[()]"))
  round_trip$replica <- as.numeric(round_trip$replica)
  
  rt_report <- round_trip |> 
    mutate(replica=as.numeric(replica)) |> 
    group_by(id,alg,replica) |> 
    summarise(avg.rt=mean(round_trips,na.rm=T)) |> ungroup() |> 
    arrange(id,replica)
  
  #Horizontal
  horizontal_report_rt <- rt_report |> select(-id) |> 
    pivot_wider(names_from = replica,values_from = avg.rt)
  # View(horizontal_report_rt)
  
  #Vertical
  vertical_report_rt <- rt_report |> select(-id) |> 
    pivot_wider(names_from = alg,values_from = avg.rt)
  # View(vertical_report_rt)
  
  table_rt_report <- vertical_report_rt
  # table_rt_report <- horizontal_report_rt
  
  boxplot_rt <- round_trip |> 
    ggplot(aes(x=factor(replica),y=round_trips))+
    geom_boxplot()+
    facet_wrap(~alg)
  boxplot_rt
  
  #Export table with averages
  jpeg(file.path(export_path,paste0(export_file_name,"_round_trips",".jpg")),width=120*ncol(table_rt_report),height =30*nrow(table_rt_report),pointsize = 30)
  grid.arrange(tableGrob(table_rt_report))
  dev.off()  
  
  # Export boxplot
  jpeg(file.path(export_path,paste0(export_file_name,"_rt_boxplot",".jpg")),width=1200,height =600,pointsize = 30)
  print(boxplot_rt)
  dev.off()  
  
}


##### REPORT: Max bounds ##### 
{
  log_bounds$id <- as.numeric(str_extract(log_bounds$alg, "\\(\\d+\\)") |> 
                                str_remove_all("[()]"))
  log_bounds$replica <- as.numeric(log_bounds$replica)
  
  bound_report <- log_bounds |> 
    mutate(replica=as.numeric(replica)) |> 
    group_by(id,alg,replica) |> 
    summarise(avg.bound=mean(log_bound,na.rm=T)) |> ungroup() |> 
    arrange(id,replica)
  
  #Horizontal
  horizontal_report_bound <- bound_report |> filter(grepl("sq",alg)) |> 
    select(-id) |> 
    pivot_wider(names_from = replica,values_from = avg.bound)
  
  
  #Vertical
  vertical_report_bound <- bound_report |> filter(grepl("sq",alg)) |> 
    select(-id) |> 
    pivot_wider(names_from = alg,values_from = avg.bound)
  
  
  table_report_bound <- vertical_report_bound
  # table_rt_report <- horizontal_report_bound
  
  boxplot_bounds <- log_bounds |> 
    filter(grepl("sq",alg)) |> 
    ggplot(aes(x=factor(replica),y=log_bound))+
    geom_boxplot()+
    facet_wrap(~alg)
  boxplot_bounds
  
  #Export table with averages
  jpeg(file.path(export_path,paste0(export_file_name,"_log_bounds",".jpg")),width=120*ncol(table_report_bound),height =30*nrow(table_report_bound),pointsize = 30)
  grid.arrange(tableGrob(table_report_bound))
  dev.off()  
  
  # Export boxplot
  jpeg(file.path(export_path,paste0(export_file_name,"_bound_boxplot",".jpg")),width=1200,height =600,pointsize = 30)
  print(boxplot_bounds)
  dev.off()  
}
  
  
  
  
  
  
  
  
  
  
  
  




