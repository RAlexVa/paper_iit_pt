rm(list=ls())
# Import libraries
library(stringr)
library(tidyverse)
library(gridExtra)# For tables
library(latex2exp) #For using latex
library(survminer); library(survival) #For plot showing how fast each replica visits all modes


# current_theme <- theme_get()
theme_set(theme_minimal());
theme_update(base_size = 17,legend.key.size = unit(1, 'cm'));

### Process simulation results ###

# Choose dimension
chosen_dim <- "highdim"; file_dim <- "highd"
# chosen_dim <- "lowdim";file_dim <- "lowd" #,10000,1000,5000
print_bimodal <- FALSE
print_multimodal <- FALSE
# chosen_ids <-c(730,733,736,739)+1
chosen_ids <-c(745,748,751,759)#+1 #For results with modes separated (initial setup, dim 1k)
chosen_ids <- c(767,769:773)#c(766,768)#c(767,769:773) #For results with modes overlaping
chosen_ids <- c(782,783,784,785)#For spaced problem
chosen_ids <- c(806:809)#For MultiM dim 5k theta 0.001
chosen_ids <- c(802:805)#For MultiM dim 3k theta 0.001
chosen_ids <- c(1024:1035)#For MultiM dim 1k varios theta
chosen_ids <- c(858:861,866:869,874:877)
chosen_ids <- c(802,804,878,879)#Para comparar lo que ya habia salido con este
chosen_ids <- c(806,808,862,864,880:881)#Para comparar lo que ya habia salido con este en dim 5k

# chosen_ids <-755
# chosen_ids <- c(650:661,663:669)
#### Chosen for lowdim bimodal problem ####
#We stay with 3 temperatures for everything
#For PT-IIT the last temperature does not achieve the 0.23 swap rate (it gets bigger rate)
#For PT A-IITm there doesn't seem to be a big issue
#For PT A-IITw it seems to have been better with 4 temperatures
# but still doesn't get the needed swap rate
# 
# chosen_bimodal <- c(202,203,282,252,254)#c(202,203,204,233,235,272,282,252,254)#c(129,135,137,139)
# print_bimodal <- TRUE
# chosen_ids <- chosen_bimodal

#### Chosen for lowdim multimodal problem ####
# For PT-IIT and PT A-IITw we only use 3 temperatures because it had the best performance in TVD
# For PT A-IITm we still have to identify the best temperature
# 
# chosen_multimodal <- c(205,206,328,294)#c(205,206,207,237,239,319,328,294)#c(165,167,169,190)
# print_multimodal <- TRUE
# chosen_ids <- chosen_multimodal


#List of files
parameters <- read_csv(paste0("inputs/simulation_details_",file_dim,".csv"), col_types = cols())
#Create table with available files
# data_sum <- tibble(file_names=list.files(path = "results", pattern = "^sim_.*\\.Rds")) |> 
#   mutate(id=as.numeric(str_extract(file_names, "(?<=id_)[0-9]+(?=\\.Rds)")),
#          dim=str_extract(file_names, "(?<=sim_)[^_]+(?=_id)")) |> 
#   filter(dim==chosen_dim) |> 
#   left_join(parameters, by="id")
  
data_sum <- tibble(file_names=list.files(path = "results", pattern = "^sim_.*\\.Rds")) |> 
  mutate(id=as.numeric(str_extract(file_names, "id_([0-9]+)",group=1)),
         sim_id=str_extract(file_names, "id_[0-9]+_([0-9]+)\\.Rds",group=1),
         dim=str_extract(file_names, "(?<=sim_)[^_]+(?=_id)")) |> 
  filter(dim==chosen_dim) |> 
  left_join(parameters, by="id")



# filter IDs to compare
data_sum <- data_sum |> filter(id %in% chosen_ids)

if(chosen_dim=="lowdim"){
Rcpp::sourceCpp("functions/cpp_functions.cpp") #To use vec_to_num function
source("functions/r_functions.R")  
check_number_modes <- unique(data_sum |> pull(model))
if(length(check_number_modes)==1){only_1_model <- TRUE;}else{only_1_model <- FALSE}
if(!only_1_model){stop("You have low_dim models with different number of modes")}else{
# Check number of replicas
  check_number_replicas <- data_sum |> select(id,matches("^t\\d+")) #Select all temperatures columns
  
  #Check which columns are full of NAs to exclude them
  na_temps <- data_sum |> 
    select(matches("^t\\d+$")) |>  # Select columns of the form tX
    summarize(across(everything(), ~ all(is.na(.x)))) |>   # Check if all values are NA
    pivot_longer(everything(), names_to = "column", values_to = "is_all_na") |>   # Convert to long format
    filter(is_all_na) |>   # Keep only columns where is_all_na is TRUE
    pull(column) 
  
  na_bf <- paste0("bf",na_temps |> str_extract("\\d+"))
  check_number_replicas <- check_number_replicas|> select(-all_of(na_temps))
  check_dif_num_replica <- any(is.na(check_number_replicas))
  
  
  if(check_dif_num_replica){#In case there's some ids with different number of replicas
    ind_rep <- which(is.na(check_number_replicas),T)[1]
    writeLines(paste0("id: ",check_number_replicas[ind_rep,1]," has less replicas"))
    stop("Some of the chosen simulations have different number of replicas")
  }else{
    num_replicas <- ncol(check_number_replicas |> select(-id))
  }
  
  
  # to store time to visit modes
  mode_time <- as.data.frame(matrix(nrow=0,ncol=4)); colnames(mode_time) <- c("alg","sim","mode","time")
  # to store time to visit modes
  tvd_report <- as.data.frame(matrix(nrow=0,ncol=5));colnames(tvd_report) <- c("alg","sim","measurement","time","tvd")
  if(check_number_modes=="7_mode"){
    #Low dim is the example with 7 modes and 4 temperatures
    tvd <- data.frame(alg=character(),sim=numeric(),tvd=numeric())
    mode_visit <- as.data.frame(matrix(nrow=0,ncol=10)); colnames(mode_visit) <- c("alg","sim","interswap",1:7)
    round_trip <- as.data.frame(matrix(nrow=0,ncol=(num_replicas+2))); colnames(round_trip) <- c("alg","sim",1:num_replicas)
    swap_rate <- as.data.frame(matrix(nrow=0,ncol=(num_replicas+1))); colnames(swap_rate) <- c("alg","sim",1:(num_replicas-1))
    iterations <- as.data.frame(matrix(nrow=0,ncol=(num_replicas+2))); colnames(iterations) <- c("alg","sim",1:num_replicas) #For the 4 replicas
    pi_modes <- as.data.frame(matrix(nrow=0,ncol=9));colnames(pi_modes) <- c("alg","sim",1:7)
    # Low dimensional true probability setup
    {
      Rcpp::sourceCpp("functions/cpp_functions.cpp") #To use vec_to_num function
      source("functions/r_functions.R")
      p <- 16 #dimension
      theta <- 15 #tail weight parameter
      
      # Modes definition
      mod1 <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)#rep(1,p)
      mod2 <- c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0)
      mod3 <- c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1)
      mod4 <- c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0)
      mod5 <- c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1)
      mod6 <- c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)
      mod7 <- c(0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0)
      modes <- c(vec_to_num(mod1),  vec_to_num(mod2),  vec_to_num(mod3),  vec_to_num(mod4),  vec_to_num(mod5),  vec_to_num(mod6),  vec_to_num(mod7))
      modes_list <- list(mod1,mod2,mod3,mod4,mod5,mod6,mod7)
      
      ### Compute true probability
      pi.true <- rep(0,2^p)
      for(i in 0:(2^p -1)){
        pi.true[i+1] <-  lik_comp(NumberToVector(i,p),modes_list,theta,p)
      }
      pi.true <- pi.true/(length(modes_list)*(1+exp(-theta))^p)
    }
    
  }
  if(check_number_modes=="bimodal"){
    #Low dim is the example with 2 modes and 4 temperatures
    tvd <- data.frame(alg=character(),sim=numeric(),tvd=numeric())
    mode_visit <- as.data.frame(matrix(nrow=0,ncol=5)); colnames(mode_visit) <- c("alg","sim","interswap",1:2)
    round_trip <- as.data.frame(matrix(nrow=0,ncol=(num_replicas+2))); colnames(round_trip) <- c("alg","sim",1:num_replicas)
    swap_rate <- as.data.frame(matrix(nrow=0,ncol=(num_replicas+1))); colnames(swap_rate) <- c("alg","sim",1:(num_replicas-1))
    iterations <- as.data.frame(matrix(nrow=0,ncol=(num_replicas+2))); colnames(iterations) <- c("alg","sim",1:num_replicas) #For the 4 replicas
    pi_modes <- as.data.frame(matrix(nrow=0,ncol=4));colnames(pi_modes) <- c("alg","sim",1:2)
    
    ##### Low-dimensional multimodal setup #####
    {
      Rcpp::sourceCpp("functions/cpp_functions.cpp") #To use vec_to_num function
      source("functions/r_functions.R")
      p <- 16 #dimension
      theta <- c()
      theta[1] <- 6 #tail weight parameter
      theta[2] <- 6 #tail weight parameter for second mode
      # Modes definition
      mod2 <- c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0)
      mod3 <- c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1)
      modes <- c(vec_to_num(mod2),  vec_to_num(mod3))
      modes_list <- list(mod2,mod3)
      
      ### Compute true probability
      pi.true <- rep(0,2^p)
      for(i in 0:(2^p -1)){
        pi.true[i+1] <-  lik_comp(NumberToVector(i,p),modes_list,theta,p)
      }
      pi.true <- pi.true/((1+exp(-theta[1]))^p + (1+exp(-theta[2]))^p)
      sum(head(sort(-pi.true),n=20))
      pi.true[modes+1]
      sum(pi.true[modes+1])
    }
  }
}
}
if(chosen_dim=="highdim"){
  Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp") #To use eval_lik function
  source("functions/r_functions.R") #To get dimension of the problem
  #High dim we can have multiple temperatures
  #we dont track distribution convergence
  #we track visit to high probability states
  
  check_number_replicas <- data_sum |> select(id,matches("^t\\d+")) #Select all temperatures columns
  # check_number_sim <- unique(data_sum |> pull(simulations))
  check_number_sim <- 2
  if(length(check_number_sim)>1){stop("Some IDs have different number of simualations")}
  #Check which columns are full of NAs to exclude them
  na_temps <- data_sum |> 
    select(matches("^t\\d+$")) |>  # Select columns of the form tX
    summarize(across(everything(), ~ all(is.na(.x)))) |>   # Check if all values are NA
    pivot_longer(everything(), names_to = "column", values_to = "is_all_na") |>   # Convert to long format
    filter(is_all_na) |>   # Keep only columns where is_all_na is TRUE
    pull(column) 
  
  na_bf <- paste0("bf",na_temps |> str_extract("\\d+"))
  check_number_replicas <- check_number_replicas|> select(-all_of(na_temps))
  check_dif_num_replica <- any(is.na(check_number_replicas))
  
  
  if(check_dif_num_replica){#In case there's some ids with different number of replicas
    ind_rep <- which(is.na(check_number_replicas),T)[1]
    writeLines(paste0("id: ",check_number_replicas[ind_rep,1]," has less replicas"))
    stop("Some of the chosen simulations have different number of replicas")
  }else{
    num_replicas <- ncol(check_number_replicas |> select(-id))
    round_trip <- as.data.frame(matrix(ncol=(num_replicas+2),nrow=0)); colnames(round_trip) <- c("alg","sim",1:num_replicas)
    if(num_replicas>1){swap_rate <- as.data.frame(matrix(ncol=(num_replicas+1),nrow=0)); colnames(swap_rate) <- c("alg","sim",1:(num_replicas-1))}
    iterations <- as.data.frame(matrix(ncol=(num_replicas+2),nrow=0)); colnames(iterations) <- c("alg","sim",1:num_replicas)
    
    list_of_states <- list()
    iter_visit <- as.data.frame(matrix(ncol=max(data_sum$simulations)+2,nrow=0));colnames(iter_visit) <- c("alg","state",1:max(data_sum$simulations))
    loglik_visited <- as.data.frame(matrix(nrow=0,ncol=4));colnames(loglik_visited) <- c("alg","sim","state","loglik")
    max_lik <- as.data.frame(matrix(nrow=0,ncol=4));colnames(max_lik) <- c("alg","sim","iteration","lik")
    
    if(check_number_sim==1){
      distances <- tibble(alg=character(),sim=numeric(),temperature=numeric(),iteration=numeric(),mode=character())
    }else{
      distances <- tibble(alg=character(),sim=numeric(),temperature=numeric(),mode=character(),min_dist=numeric(),max_dist=numeric(),min_iter=numeric(),max_iter=numeric())  
    }
  }
}
time_taken <- as.data.frame(matrix(nrow=0,ncol=2));colnames(time_taken) <- c("alg","time")
full_iter_names <- c() #To identify the list of full_iter with the corresponding algorithms
full_iter <- list()
k <- 1
Q <- 1

# Start creating datasets with information
for(i in 1:nrow(data_sum)){
  data <- readRDS(file.path("results",data_sum[i,1]))
  selected_id <- data_sum |> slice(i)|> pull(id)#Extract ID
  tot_sim <- data_sum |> slice(i)|> pull(simulations)
  algorithm <- data_sum |> slice(i) |> pull(algorithm)
  tot_iter <- data_sum |> slice(i) |> pull(iterations)
  tot_swap <- data_sum |> slice(i) |> pull(total_swap)
  interswap <- data_sum |> slice(i) |> pull(interswap)
  temperatures <- as.numeric(data_sum |> slice(i) |> select(matches("^t\\d{1,2}$")))
  temperatures <- temperatures[!is.na(temperatures)]# all the temperatures are in order and consecutive in the CSV file
  num_modes <- data_sum |> slice(i) |> pull(num_modes)
  if(algorithm=="PT_A_IIT"){algorithm <- paste0("PT A-IIT (",data_sum |> slice(i)|> pull(bf),")")}
  if(algorithm=="PT_IIT_no_Z"){algorithm <- "PT-IIT no Z"}
  if(algorithm=="PT_IIT_Z"){algorithm <- paste0("PT-IIT (",data_sum |> slice(i)|> pull(bf),")")}
  if(algorithm=="PT_A_IIT_RF"){algorithm <- "PT A-IIT w"}
  if(chosen_dim=="lowdim"){  if(data_sum |> slice(i) |> pull(reduc_model)=="zero"){algorithm <- "RF-MH"}}
  
### Option for when the Number of simulations is in the name
  simid_in_name <- !is.na(data_sum |> slice(i) |> pull(sim_id))
  if(simid_in_name){
    list_sim_numbers <- as.numeric(data_sum |> slice(i) |> pull(sim_id))
  }else{
    list_sim_numbers <- 1:tot_sim
  }
  
##### Optional add the ID of the simulation into the name of the algorithm
  if(!(print_bimodal || print_multimodal)){algorithm <- paste0(algorithm,"(",data_sum |> slice(i) |> pull(id),")")}
  print(data_sum[i,"algorithm"])
  print(names(data))

#Extract time
temp <- as.data.frame(data[["time_taken"]])
if(!is_empty(temp)){
  temp$alg <- algorithm
  temp <- temp |> select(alg,everything())
  colnames(temp) <- c("alg","time")
  time_taken <- rbind(time_taken,temp)
}

### Specific extractions for lowdim example
if(chosen_dim=="lowdim"){
  # Extract TVD
  temp <- tibble(alg=algorithm,sim=1:tot_sim,tvd=data[["tvd"]])
  tvd <- rbind(tvd,temp)
  
  # Extract visit of modes
  temp <- as.data.frame(data[["mode_visit"]])
  colnames(temp) <- 1:ncol(temp)
  temp$sim <- list_sim_numbers
  temp$alg <- algorithm
  temp$interswap <- data_sum |> slice(i) |> pull(interswap)
  temp <- temp |> select(alg,sim,interswap,everything())
  mode_visit <- rbind(mode_visit,temp)
  
  #Extract dist. estimation for modes
  temp <- data[["est_pi"]]
  if(max(colSums(temp))>1){  #check if distribution estimate is not normalized
    temp <- temp/colSums(temp)
  }
  temp <- as.data.frame(t(temp[modes+1,]))
  colnames(temp) <- 1:ncol(temp)
  temp$sim <- list_sim_numbers
  temp$alg <- algorithm
  temp <- temp |> select(alg,sim,everything())
  pi_modes <- rbind(pi_modes,temp)

  #Extract time to visit each of the modes
  temp <- data[["time_visit"]]
  if(!is_empty(temp)){
    temp <- as.data.frame(t(temp[modes+1,])) #As many rows as simulations as many columns as modes
    #define colnames
    if(ncol(temp)==2){
      colnames(temp) <- c("mod2","mod3")
    }else   if(ncol(temp)==7){
      colnames(temp) <- paste0("mod",1:7)
    }else{stop("Number of modes in time_visit is not 7 or 2")}
    
    temp$sim <- list_sim_numbers
    temp$alg <- algorithm
    temp <- temp |> select(alg,sim,everything()) |> pivot_longer(starts_with("mod"),names_to = "mode",values_to="time")
    mode_time <- rbind(mode_time,temp)
  }else{
    print(paste0("ID: ",data_sum[i,"id"],"doesn't have time_visit"))
  }
  #Extract time to TVD measurement
  temp <- data[["tvd_time_report"]]
  if(!is_empty(temp)){
    # colnames(temp) <- paste0("t",1:ncol(temp))
    temp <- as.data.frame(temp)
    temp$sim <- list_sim_numbers
    temp$alg <- algorithm
    temp <- temp |> select(alg,sim,everything())
    temp <- temp |> pivot_longer(cols=starts_with("V"),names_to="measurement",values_to="time")
    
    temp2 <- data[["tvd_report"]]
    temp2 <- as.data.frame(temp2)
    temp2$sim <- list_sim_numbers
    temp2$alg <- algorithm
    temp2 <- temp2 |> select(alg,sim,everything())
    temp2 <- temp2 |> pivot_longer(cols=starts_with("V"),names_to="measurement",values_to="tvd")
    
    temp <- left_join(temp,temp2,by=c("alg","sim","measurement"))
    
    temp$measurement <- as.numeric(gsub("V","",temp$measurement))
    
    tvd_report <- rbind(tvd_report,temp)
  }
  
}
if(chosen_dim=="highdim"){
### Specific extractions for highdim example
      p <- data_sum$p;

##### Extract distance to modes
    output_name <- paste0("distance_modes")
    output_time <- paste0("time_modes")
    
    ### Extract minimum distances  
    temp_m <- as.data.frame(data[[output_name]]) 
    colnames(temp_m) <- round(temperatures,2) #column name is the temperature value
    temp_m$alg <- algorithm #Column with name of the algorithm
    temp_m$mode <- paste0("m",1:num_modes) #Column indicating the mode
    temp_m$sim <- (list_sim_numbers) #Column indicating the number of simulation
    # temp_time_m <- as.data.frame(data[[output_name]])
    #Pivot wider
    temp_m <- temp_m |> pivot_longer(-(alg:sim),names_to="temperature",values_to = "min_dist")
    
    ### Extract time to reach minimum distances
    temporal_time <- as.data.frame(data[[output_time]])
    colnames(temporal_time) <- round(temperatures,2)
    #Create columns to identify
    temporal_time$alg <- algorithm
    temporal_time$mode <- paste0("m",1:num_modes)
    temporal_time$sim <- (list_sim_numbers)
    #Pivot longer
    temporal_time <- temporal_time |> pivot_longer(-(alg:sim),names_to="temperature",values_to = "time_find")
    
    temp_join <- left_join(temp_m,temporal_time,by=c("alg","mode","sim","temperature"))
    
    distances <- rbind(distances,temp_join)
  
  }
  
  
  if(!is_empty(data[["round_trips"]]) && ncol(data[["round_trips"]])>1){
    #Extract number of round trips rate
    temp <- as.data.frame(data[["round_trips"]])
    colnames(temp) <- 1:ncol(temp)
    temp$sim <- list_sim_numbers
    temp$alg <- algorithm
    temp <- temp |> select(alg,sim,everything())
    round_trip <- rbind(round_trip,temp)
    # Extract replica swap rate
    temp <- as.data.frame(data[["swap_rate"]])
    colnames(temp) <- 1:ncol(temp)
    temp$sim <- list_sim_numbers
    temp$alg <- algorithm
    temp <- temp |> select(alg,sim,everything())
    swap_rate <- rbind(swap_rate,temp)
  }
  if(!is_empty(data[["total_iter"]])&& ncol(data[["total_iter"]])>1){ 
    # Extract total iterations
    # dims<- dim(data[["total_iter"]])
    # full_iter[[k]] <- data[["total_iter"]]
    # full_iter_names[k] <- algorithm
    # k <- k+1;
    # temp <- as.data.frame(t(colSums(data[["total_iter"]])))/dims[1]
    # #temp is the average number of Rejection Free steps before trying a swap
    # colnames(temp) <- 1:ncol(temp)
    # temp$sim <- list_sim_numbers
    # temp$alg <- algorithm
    # temp <- temp |> select(alg,sim,everything())
    # iterations <- rbind(iterations,temp)

### Considering that we only have 1 simulation per read dataset
    temp_single <- data[["total_iter"]][,,1]
    check_rows <- rowSums(temp_single)
    cutoff <- min(which(check_rows==0))
    cutoff <- min(cutoff,length(check_rows))#In case we don't have any zeroes = didn't break the for loop
    temp_single <- temp_single[1:(cutoff-1),]
    temp <- data.frame(matrix(colSums(temp_single)/nrow(temp_single),nrow=1))
    colnames(temp) <- 1:length(temp)
    temp$sim <- list_sim_numbers
    temp$alg <- algorithm
    temp <- temp |> select(alg,sim,everything())
    iterations <- rbind(iterations,temp)
      
  }
}


##### Export plots and tables #####
export_path <- paste0("C:/Users/ralex/Documents/src/paper_iit_pt/images/",chosen_dim,"_ex")
export_file_name <- paste0(paste0(chosen_ids,collapse="_"),"_",chosen_dim)
if(print_bimodal){export_file_name <- "bimodal"}
if(print_multimodal){export_file_name <- "multimodal"}
# full_path <- file.path(export_path,paste0("tvd_",export_file_name,".jpg"))

#################################### REPORTS FOR LOW AND HIGH DIMENSION ######

##### Time taken #####

# time_plot <- time_taken|>  filter(alg!='IIT') |> 
#   ggplot(aes(x=alg,y=time,fill=alg)) +
#   geom_boxplot(show.legend = FALSE)+
#   labs(fill='Algortihm',x="",y="Time in seconds")+
#   theme_minimal(base_size = 17)+
#   theme(legend.key.size = unit(1, 'cm'))
# time_plot
time_table <- time_taken |> filter(!str_starts(alg,'IIT')) |> 
  group_by(alg) |> 
  summarise(min=min(time),
            q1=quantile(time,probs=0.25),
            median=quantile(time,probs=0.5),
            mean=mean(time),
            q3=quantile(time,probs=0.75),
            max=max(time))

jpeg(file.path(export_path,paste0(export_file_name,"_time_table",".jpg")),width=80*ncol(time_table),height=35*nrow(time_table),pointsize = 30)
grid.arrange(tableGrob(time_table))
dev.off()

(time_boxplot <- time_taken |> 
  filter(!str_starts(alg,'IIT')) |>
  ggplot(aes(x=alg,y=time))+
  geom_boxplot())

##### Report on average swap rate

#First create column names depending on the number of replicas
if(num_replicas>1){
  
  swap_names <- c()
  for(i in 1:(num_replicas-1)){
    swap_names <- c(swap_names,paste0(i,"↔",i+1))
  }
  
  swap_report <- swap_rate |> 
    group_by(alg) |>
    summarize(across(-sim, function(x) mean(x, na.rm = TRUE)))
  colnames(swap_report) <- c("alg",swap_names)
  
  jpeg(file.path(export_path,paste0(export_file_name,"_table_swap_rate",".jpg")),width=93*ncol(swap_report),height=33*nrow(swap_report),pointsize = 30)
  grid.arrange(tableGrob(swap_report))
  dev.off()
  
  ##### Report on average round trip rate
  rt_report <- round_trip |> 
    group_by(alg) |>
    summarize(across(-sim, function(x) mean(x, na.rm = TRUE)))
  colnames(rt_report) <- c("alg",paste0("R",1:num_replicas))
  
  jpeg(file.path(export_path,paste0(export_file_name,"_table_roundtrip_rate",".jpg")),width=90*ncol(rt_report),height=40*nrow(rt_report),pointsize = 30)
  grid.arrange(tableGrob(rt_report))
  dev.off()
  
}


if(chosen_dim=="lowdim"){
  ##### Delete first row with NA#####
  tvd <-  tvd |> filter(!is.na(alg)) 
  mode_visit <- mode_visit |> filter(!is.na(alg))
  pi_modes <- pi_modes|> filter(!is.na(alg))
  mode_time <- mode_time |> filter(!is.na(alg))
##### Compare estimation of modes #####  
  
  pi_modes |> pivot_longer(cols = -(alg:sim), names_to = "mode", values_to = "pi_est") |> 
    ggplot(aes(x=mode,y=pi_est,fill=alg))+
    geom_boxplot(show.legend = FALSE)+
    geom_hline(yintercept = pi.true[modes[1]+1], color = "red", linetype = "dashed", size = 1)+
    facet_wrap(~alg)
  # +
  #   theme_minimal(base_size = 17)+
  #   theme(legend.key.size = unit(1, 'cm'))
  # 
#### TVD computed only with the modes
  col_selected <- colnames(pi_modes)
  col_selected <- col_selected[col_selected!="alg" & col_selected!="sim"]
  pi_modes |> rowwise() |> 
    mutate(tvd=0.5*sum(abs(pi.true[modes+1]-c_across(col_selected)))) |> 
    ungroup() |> 
    select(alg,tvd) |> 
    ggplot(aes(x=alg,y=tvd,fill=alg)) +
    geom_boxplot(show.legend = FALSE)+
    labs(fill='Algortihm',x="",y="Total Variation Distance")
  # +
  #   theme_minimal(base_size = 17)+
  #   theme(legend.key.size = unit(1, 'cm'))

#Min and max values in TVD
 min_tvd_1 <-  pi_modes |> rowwise() |> 
    mutate(tvd=0.5*sum(abs(pi.true[modes+1]-c_across(col_selected)))) |> 
    select(alg,tvd) |> 
    group_by(alg) |> slice(which.min(tvd))
  
 max_tvd_1 <-  pi_modes |> rowwise() |> 
    mutate(tvd=0.5*sum(abs(pi.true[modes+1]-c_across(col_selected)))) |> 
    select(alg,tvd) |> 
    group_by(alg) |> slice(which.max(tvd))
  
  
min_tvd_2 <- tvd |> select(alg,tvd) |> group_by(alg) |> slice(which.min(tvd))
max_tvd_2 <- tvd |> select(alg,tvd) |> group_by(alg) |> slice(which.max(tvd))
  
min_tvd_compare <- left_join(min_tvd_1,min_tvd_2,by=c("alg"),suffix = c("_modes","_full"))
max_tvd_compare <- left_join(max_tvd_1,max_tvd_2,by=c("alg"),suffix = c("_modes","_full"))
grid.arrange(tableGrob(min_tvd_compare))
grid.arrange(tableGrob(max_tvd_compare))
##### Total Variation Distance #####
  
  tvd_plot <- tvd |>  filter(!str_starts(alg,'IIT')) |> 
    ggplot(aes(x=alg,y=tvd,fill=alg)) +
    geom_boxplot()+
    labs(fill='Algortihm',x="",y="Total Variation Distance")
# +
#     theme_minimal(base_size = 17)+
#     theme(legend.key.size = unit(1, 'cm'))
  tvd_plot
  
  jpeg(file.path(export_path,paste0(export_file_name,"_tvd",".jpg")),width=800,height =400,pointsize = 30)
  print(tvd_plot)
  dev.off()
##### First visit to modes #####
col_names <- col_selected

mode_sum <- mode_visit |> 
  rowwise()|> 
  mutate(last_visit=max(c_across(all_of(col_names))), first_visit=min(c_across(all_of(col_names))[c_across(all_of(col_names))>0])) |> 
  mutate(first_mode =  names(pick(all_of(col_names)))[which(c_across(all_of(col_names)) == first_visit)[1]]) |> 
  mutate(last_mode =  names(pick(all_of(col_names)))[which(c_across(all_of(col_names)) == last_visit)[1]]) |> 
  mutate(total_modes = sum(c_across(all_of(col_names)) > 0)) |> 
  select(-all_of(col_names))

##### Report on number of modes visited by each algorithm 
table_visited <- mode_sum |> rename(algorithm=alg) |> 
  group_by(algorithm,total_modes) |> 
  summarise(count=n()) |> ungroup() |> arrange(total_modes,desc(algorithm)) |> 
  pivot_wider(id_cols =algorithm,
              names_from = total_modes,
              values_from = count,
              values_fill=0)

jpeg(file.path(export_path,paste0(export_file_name,"_table_visited_modes",".jpg")),width=100 + (60*(ncol(table_visited)-1)),height=40*nrow(table_visited),pointsize = 30)
grid.arrange(tableGrob(table_visited))
dev.off()

##### Report on number of iterations for the original replica to visit all modes (most of the times after a swap)
iterations_to_explore <- mode_sum |> filter(!str_starts(alg,'IIT')) |> 
  group_by(alg) |> 
  summarise(min=min(last_visit),
            q1=quantile(last_visit,probs=0.25),
            median=quantile(last_visit,probs=0.5),
            mean=mean(last_visit),
            q3=quantile(last_visit,probs=0.75),
            max=max(last_visit))

jpeg(file.path(export_path,paste0(export_file_name,"_table_iterations",".jpg")),width=60*ncol(iterations_to_explore),height=35*nrow(iterations_to_explore),pointsize = 30)
grid.arrange(tableGrob(iterations_to_explore))
dev.off()


##### Report on number of replica swaps needed to visit all modes
# First algorithms do 2k iterations before trying a replica swap
swaps_to_explore <- mode_sum |> filter(!str_starts(alg,'IIT'),!str_starts(alg,'PT A-IIT m')) |> 
  mutate(last_visit=last_visit/interswap) |> 
  group_by(alg) |> 
  summarise(min=min(last_visit),
            q1=quantile(last_visit,probs=0.25),
            median=quantile(last_visit,probs=0.5),
            mean=mean(last_visit),
            q3=quantile(last_visit,probs=0.75),
            max=max(last_visit))

if(nrow(iterations)>0){
  temp <- mode_sum |> filter(str_starts(alg,'PT A-IIT m')) |>
    select(alg,sim,last_visit) |> 
    left_join(iterations,by=c("alg","sim")) |> 
    mutate(last_visit=last_visit/`1`) |> 
    select(alg,last_visit) |> 
    group_by(alg) |> 
    summarise(min=round(min(last_visit),2),
              q1=round(quantile(last_visit,probs=0.25),2),
              median=round(quantile(last_visit,probs=0.5),2),
              mean=round(mean(last_visit),2),
              q3=round(quantile(last_visit,probs=0.75),2),
              max=round(max(last_visit),2))
  swaps_to_explore <- rbind(swaps_to_explore,temp)
  
## Report on number of iterations interswap
  
  interswap_report <- data_sum |> select(id,interswap)
  
  rep_iter <- iterations |> 
    mutate(id=as.numeric(str_extract(alg,"\\d+"))) |> 
    left_join(interswap_report,by="id") |> 
    select(-id) |> 
    group_by(alg,interswap) |> 
    summarise(across(-sim, mean))

  if(print_bimodal || print_multimodal){rep_iter <- rep_iter |> select(-interswap)}
  jpeg(file.path(export_path,paste0(export_file_name,"_interswaps",".jpg")),width=85*ncol(rep_iter),height=35*nrow(rep_iter),pointsize = 30)
  grid.arrange(tableGrob(rep_iter))
  dev.off()  
  
  
}

jpeg(file.path(export_path,paste0(export_file_name,"_table_swaps",".jpg")),width=100*ncol(swaps_to_explore),height=35*nrow(swaps_to_explore),pointsize = 30)
grid.arrange(tableGrob(swaps_to_explore))
dev.off()

#### Report on time to visit the modes  

(plot_time_mode <- mode_time |> group_by(alg,sim) |>
  summarise(first_time=min(time),last_time=max(time)) |> 
  ggplot(aes(x=alg,y=last_time, fill=alg)) +
  geom_boxplot()+
  labs(fill='Algortihm',y="seconds",x="",title="Time to find last mode")+
    theme(axis.text.x = element_blank()))
  # theme_minimal(base_size = 17)+
  # theme(legend.key.size = unit(1, 'cm'),
  #       axis.text.x = element_blank()))

mode_time |> group_by(alg,sim) |>
  summarise(first_time=min(time),last_time=max(time)) |> ungroup() |> 
  ggplot(aes(x=first_time,y=last_time, col=alg)) +
  geom_point()

forsurv <- mode_time |>group_by(alg,sim) |>
  summarise(first_time=min(time),last_time=max(time)) |> ungroup() |> 
  select(alg,last_time)

# forsurv <- forsurv |> filter(str_detect(alg,"\\(26[7]\\)|\\(27[1268]\\)|\\(28[023]\\)"))

fit <- survfit(Surv(last_time,rep(1,nrow(forsurv)))~alg,data=forsurv)
if(check_number_modes=="bimodal"){time_br <- 0.2}
if(check_number_modes=="7_mode"){time_br <- 0.5}
(plot_surv_mode <- ggsurvplot(fit,
           data=forsurv,
           fun="event",
           palette = "Set1",    # Color palette
           xlab = "Time (seconds)",
           ylab = "Prop. of simulations visiting all modes",
           legend.title = "Algorithm",
           # break.time.by = time_br,
           font.x = 15,        # X-axis label font size
           font.y = 15,        # Y-axis label font size
           font.tickslab = 12, # Axis tick labels (numbers) font size
           font.legend = 10))   # Legend text font size)

table_time_mode <- mode_time |> group_by(alg,sim) |>
  summarise(first_time=min(time),last_time=max(time)) |> 
  ungroup() |> group_by(alg) |> 
summarise(min=min(last_time),
          q1=quantile(last_time,probs=0.25),
          median=quantile(last_time,probs=0.5),
          mean=mean(last_time),
          q3=quantile(last_time,probs=0.75),
          max=max(last_time))

grid.arrange(tableGrob(table_time_mode))

jpeg(file.path(export_path,paste0(export_file_name,"_time_mode",".jpg")),width=800,height =400,pointsize = 30, quality=100)
print(plot_surv_mode)
dev.off()

### TVD report

sum_tvd <- tvd_report |> mutate(measurement=measurement/max(measurement)) |> 
  group_by(alg,measurement) |> 
  summarise(mean_time=mean(time),
            min_time=min(time),
            q1_time=quantile(time,0.25),
            q2_time=quantile(time,0.5),
            q3_time=quantile(time,0.75),
            max_time=max(time),
            mean_tvd=mean(tvd),
            min_tvd=min(tvd),
            q1_tvd=quantile(tvd,0.25),
            q2_tvd=quantile(tvd,0.5),
            q3_tvd=quantile(tvd,0.75),
            max_tvd=max(tvd)) 


max_tvd_plot <- sum_tvd |> 
  ggplot(aes(x=max_time,y=max_tvd, col=alg))+
  geom_point()+
  geom_line()+
  scale_x_continuous()+
  labs(title='Maximum TVD',color='Algorithm', x="Seconds", y="Total Variation Distance")
max_tvd_plot
mean_tvd_plot <- sum_tvd |> 
  ggplot(aes(x=max_time,y=mean_tvd, col=alg))+
  geom_point()+
  geom_line()+
  scale_x_continuous()+
  labs(title='Mean TVD',color='Algorithm', x="Seconds", y="Total Variation Distance")
mean_tvd_plot 

trunc_mean_tvd_plot <- sum_tvd |> filter(mean_tvd<0.1) |> 
  ggplot(aes(x=max_time,y=mean_tvd, col=alg))+
  geom_point(size=2.5)+
  geom_line(linewidth=0.8)+
  scale_x_continuous()+
  labs(title='Mean TVD',color='Algorithm', x="Seconds", y="Total Variation Distance")+
  theme_minimal()
trunc_mean_tvd_plot

sum_tvd |> 
  ggplot(aes(x=max_time,y=q2_tvd, fill=alg, col=alg))+
  geom_line(linewidth=0.8)+
  geom_point(pch=1,size=1.5)+
  geom_ribbon(aes(ymin=q1_tvd, ymax=q3_tvd),alpha=0.04)
  



sum_tvd |> 
  ggplot(aes(x=max_time,y=q2_tvd, col=alg))+
  geom_point()+
  geom_line()+
  scale_x_continuous()+
  labs(title='Median TVD',color='Algorithm', x="Seconds", y="Total Variation Distance")

sum_tvd |> 
  ggplot(aes(x=max_time,y=min_tvd, col=alg))+
  geom_point()+
  geom_line()+
  scale_x_continuous()+
  labs(title='Minimum TVD',color='Algorithm', x="Seconds", y="Total Variation Distance")



jpeg(file.path(export_path,paste0(export_file_name,"_tvd_mean",".jpg")),width=800,height =400,pointsize = 30)
print(mean_tvd_plot)
dev.off()

jpeg(file.path(export_path,paste0(export_file_name,"_tvd_mean_trunc",".jpg")),width=800,height =400,pointsize = 30)
print(trunc_mean_tvd_plot)
dev.off()

jpeg(file.path(export_path,paste0(export_file_name,"_tvd_max",".jpg")),width=800,height =400,pointsize = 30)
print(max_tvd_plot)
dev.off()


}# Finish low dim reports

#Starts high dim reports
if(chosen_dim=="highdim"){

#Input NAs
  # distances$time_find[distances$time_find==0] <- NA

## Get out ID
    distances$id <- str_extract(distances$alg, "\\(\\d+\\)") |> 
      str_remove_all("[()]")  # Extract numbers between parentheses and remove the parentheses
    distances$algorithm <- str_replace(distances$alg, "\\s*\\(.*", "")  # Extract text before parenthesis
### Speed to modes for maximum temperature
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
    # View(distances)
    View(distances |> 
           group_by(alg,mode, sim) |> 
           summarise(minim_dist=min(min_dist),max_dist=max(min_dist)) |> 
           group_by(alg,mode) |> 
           summarise(min.avg=mean(minim_dist), max.avg=mean(max_dist)))

    
    
#### Reports   
### This is for specific Ids
    
    new_id_join <- tibble(alg=as.character(c("PT A-IIT m(745)",
                                             "PT-IIT (sq)(748)",
                                             "PT-IIT (sq)(751)",
                                             "PT A-IIT m(746)",
                                             "PT-IIT (sq)(749)",
                                             "PT-IIT (sq)(752)",
                                             "PT A-IIT m(747)",
                                             "PT-IIT (sq)(750)",
                                             "PT-IIT (sq)(753)",
                                             "PT A-IIT m(759)",
                                             "PT A-IIT m(760)",
                                             "PT A-IIT m(761)",
                                             "PT A-IIT m(806)",  
                                             "PT A-IIT m(808)",  
                                             "PT-IIT (min)(809)",
                                             "PT-IIT (sq)(807)",
                                             "PT A-IIT m(802)",  
                                             "PT-IIT (sq)(803)",
                                             "PT A-IIT m(804)",  
                                             "PT-IIT (min)(805)")),
                          new_alg=c("PT_A_IIT(300)",
                                    "PT_IIT_Z(5)",
                                    "PT_IIT_Z(10)",
                                    "PT_A_IIT(300)",
                                    "PT_IIT_Z(5)",
                                    "PT_IIT_Z(10)",
                                    "PT_A_IIT(300)",
                                    "PT_IIT_Z(5)",
                                    "PT_IIT_Z(10)",
                                    "PT MH-mult (300)",
                                    "PT MH-mult (300)",
                                    "PT MH-mult (300)",
                                    "PT_A_IIT",
                                    "PT MH-mult",
                                    "PT RF-MH",
                                    "PT-IIT",
                                    "PT_A_IIT",
                                    "PT-IIT",
                                    "PT MH-mult",
                                    "PT RF-MH"))

    
#Apply modification to datasets    
    dist_t1_times <- dist_t1_times |> left_join(new_id_join,by="alg")
    swap_rate <- swap_rate|> left_join(new_id_join,by="alg")
    iterations <- iterations |> left_join(new_id_join,by="alg")
 
    dist_t1_times$new_alg <- dist_t1_times$alg
    swap_rate$new_alg <- swap_rate$alg
    iterations$new_alg <- iterations$alg
    
### Report minimum distance reached by ANY replica
report_min_dist <- distances |> 
  group_by(alg,mode) |> 
  slice_min(min_dist,n=1) |> 
  slice_min(time_find,n=1) |> 
  select(alg,mode,min_dist,time_find) |> 
  ungroup()
 view(report_min_dist)   
    
wide_report_min_dist <-  report_min_dist |> 
   select(-time_find) |> 
   pivot_wider(names_from = mode,values_from = min_dist)
 view(wide_report_min_dist)
 
#Report of simulations that ANY replica have visited the modes.

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
    
    #Add new identification of algorithm
    first_visit_report <- first_visit_report |> left_join(new_id_join,by="alg")
    ###OR
    first_visit_report$new_alg <- first_visit_report$alg
    
    
    report_modes_visited <- first_visit_report |> 
      select(new_alg,algorithm,id,sim,starts_with("min_m")) |> 
      pivot_longer(cols=starts_with("min_m"),names_to = "mode",values_to="time") |> 
      mutate(visited=time!=Inf) |> 
      group_by(id,new_alg,sim) |> 
      summarise(modes_visited=sum(visited)) |> 
      ungroup() 
    
    sim_not_visiting <- report_modes_visited |> 
      filter(modes_visited<max(report_modes_visited$modes_visited)) |> 
      pull(sim)
    
    
    unique_ids <- unique(report_modes_visited |> 
                           select(id,new_alg))
    
    sim_ids <- 101:250
    
    full_list <- tibble(id=rep(unique_ids$id,length(sim_ids)),
                        new_alg=rep(unique_ids$new_alg,length(sim_ids)),
                                    sim=rep(sim_ids,each=length(unique_ids$new_alg)))
    
    #Report on simulations that didn't finish or didin't find all modes
    number_modes <- 2
    detailed_report <- full_list |> 
      left_join(report_modes_visited,by=c("id","new_alg","sim")) |> 
      select(-new_alg) |> 
      mutate(across(-sim, ~replace_na(.,0))) |> 
      pivot_wider(names_from = id,values_from=modes_visited) |> 
      rowwise() |> 
      mutate(none_finished=all(c_across(as.character(chosen_ids))<number_modes),
             all_finished=all(c_across(as.character(chosen_ids))==number_modes),
             some_finished=any(c_across(as.character(chosen_ids))==number_modes) && !all_finished)
    
    sim_to_rerun <- detailed_report |> filter(some_finished==T,none_finished==F)
    sim_no_finish <- detailed_report |> filter(none_finished==T)
    
    jpeg(file.path(export_path,paste0(export_file_name,"sim_no_finish",".jpg")),width=110 + (60*(ncol(sim_no_finish)-1)),height=25*nrow(sim_no_finish),pointsize = 30)
    grid.arrange(tableGrob(sim_no_finish))
    dev.off()
    
    jpeg(file.path(export_path,paste0(export_file_name,"sim_to_rerun",".jpg")),width=90 + (60*(ncol(sim_to_rerun)-1)),height=30*nrow(sim_to_rerun),pointsize = 30)
    grid.arrange(tableGrob(sim_to_rerun))
    dev.off()
    
    
    jpeg(file.path(export_path,paste0(export_file_name,"_detailed_report",".jpg")),width=110 + (60*(ncol(detailed_report)-1)),height=25*nrow(detailed_report),pointsize = 30)
    grid.arrange(tableGrob(detailed_report))
    dev.off()
    
    long_report <- full_list |> left_join(report_modes_visited,by=c("id","new_alg","sim")) |> 
      filter(modes_visited<max(report_modes_visited$modes_visited) | is.na(modes_visited)) 
    
    
    jpeg(file.path(export_path,paste0(export_file_name,"_not_visiting",".jpg")),width=110 + (60*(ncol(long_report)-1)),height=25*nrow(long_report),pointsize = 30)
    grid.arrange(tableGrob(long_report))
    dev.off()
    #List of sim (seeds) that did not finish or didn't find all modes
    unique(full_report |> pull(sim))
    #Count of simulations that didnot finish or didn't find all modes
    summarized_report <- full_list |> left_join(report_modes_visited,by=c("id","new_alg","sim")) |> 
      mutate(modes_visited=ifelse(is.na(modes_visited),0,modes_visited)) |> 
      # filter(modes_visited<max(report_modes_visited$modes_visited) | is.na(modes_visited)) |> 
      group_by(id,new_alg, modes_visited) |> 
      summarise(counted=n()) |> 
      ungroup() |> 
      pivot_wider(names_from=modes_visited,values_from=counted)
    
    jpeg(file.path(export_path,paste0(export_file_name,"_number_modes_visited",".jpg")),width=50 + (100*(ncol(summarized_report)-1)),height=50*nrow(summarized_report),pointsize = 30)
    grid.arrange(tableGrob(summarized_report))
    dev.off()
    
    report_replicas_visiting <- first_visit_report |> 
      select(new_alg,algorithm,id,starts_with("row_m")) |> 
      pivot_longer(cols=starts_with("row_m"),names_to = "mode",values_to="replica") |> 
      group_by(id,new_alg,replica) |> 
      summarise(count=n()) |> 
      pivot_wider(names_from = replica,values_from = count)
    
    jpeg(file.path(export_path,paste0(export_file_name,"_replicas_visiting",".jpg")),width=110 + (60*(ncol(report_replicas_visiting)-1)),height=40*nrow(report_replicas_visiting),pointsize = 30)
    grid.arrange(tableGrob(report_replicas_visiting))
    dev.off()
    
    
    
    
    
##################################################
    #Report of speed to mode considering that any replica visits the mode
    forsurv <- first_visit_report |> select(new_alg,last_visit)
    forsurv <- forsurv |> filter(last_visit<Inf)
    # forsurv <- forsurv |> filter(new_alg!="PT_IIT_Z(5)")
    forsurv_bk <- forsurv
    ids_to_print <- chosen_ids
    ids_to_print <- c(862,864)#c(862,864)#c(854,856)#c(854:857)
    ids_to_print <- c(806,808,880,881)
    forsurv <- forsurv_bk |> filter(grepl(paste(ids_to_print,collapse = "|"),new_alg))
    fit <- survfit(Surv(last_visit,rep(1,nrow(forsurv)))~new_alg,data=forsurv)
    
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
    
    if(!identical(ids_to_print,chosen_ids)){
      export_file_name <- paste0(paste0(ids_to_print,collapse="_"),"_",chosen_dim)}
    jpeg(file.path(export_path,paste0(export_file_name,"_speed_mode_anyrep",".jpg")),width=1200,height =600,pointsize = 30)
    print(plot_surv_mode)
    dev.off()
    
    
    
############################################################          
    
    
#Report of simulations that the first temperature haven't visited all modes       
    listing_tot_sim <- dist_t1_times |> group_by(new_alg) |> summarize(tot_sim=n())
    listing_uncompleted_sim <- dist_t1_times |> filter(is.na(last_time)) |> group_by(new_alg) |> summarize(no_visit=n())
    report_uncompleted_sim <- dist_t1_times |> filter(is.na(last_time)) |> select(new_alg,sim,id) 
    report_uncompleted_sim
    
    list_ids_notcomplete <- report_uncompleted_sim  |> group_by(id) |> 
      summarize(array_string = paste(sim, collapse = ",")) |> 
      ungroup()
    list_ids_notcomplete
    report <- left_join(listing_tot_sim,listing_uncompleted_sim,by="new_alg")
    report
    
    jpeg(file.path(export_path,paste0(export_file_name,"_sim_report",".jpg")),width=110 + (60*(ncol(report)-1)),height=40*nrow(report),pointsize = 30)
    grid.arrange(tableGrob(report))
    dev.off()
# Finish report
    
    
   #Avg. swap rate report 
    swap_rate_report <- swap_rate |> select(-sim) |> 
      group_by(new_alg,alg) |> summarise(across(everything(),mean))
    
    jpeg(file.path(export_path,paste0(export_file_name,"_swap_rate_summary",".jpg")),width=110 + (80*(ncol(swap_rate_report)-1)),height=40*nrow(swap_rate_report),pointsize = 30)
    grid.arrange(tableGrob(swap_rate_report))
    dev.off()
    
    ## Average iterations interswap
    iterations_report <- iterations |> select(-sim) |> 
      group_by(new_alg,alg)|> summarise(across(everything(),mean))
    
    jpeg(file.path(export_path,paste0(export_file_name,"_iterations_summary",".jpg")),width=50 + (75*(ncol(iterations_report)-1)),height=40*nrow(iterations_report),pointsize = 30)
    grid.arrange(tableGrob(iterations_report))
    dev.off()
    
    
    boxplot_swap_rate <- swap_rate |> pivot_longer(cols=-c("alg","sim","new_alg"),names_to="replica",values_to="swap_rate") |> 
      ggplot(aes(x=replica,y=swap_rate, fill=replica))+
      geom_boxplot()+
      facet_grid(~new_alg)+ 
      theme(legend.position="none")
    
    jpeg(file.path(export_path,paste0(export_file_name,"_swap_rate_boxplot",".jpg")),width=1200,height =600,pointsize = 30)
    print(boxplot_swap_rate)
    dev.off()
   
    # new_id_join <- tibble(id=as.character(650:669),new_id=rep(c("PT_A_IIT(300)","PT_A_IIT(150)","PT_IIT_Z(100)","PT_IIT_Z(50)"),each=5))
    # forsurv <- dist_t1_times  |> left_join(new_id_join,by="id") |>
    # select(new_id,last_time)
    # fit <- survfit(Surv(last_time,rep(1,nrow(forsurv)))~new_id,data=forsurv)
 
#For bimodal with 100 simulations each
    forsurv <- dist_t1_times |> select(new_alg,last_time)
    fit <- survfit(Surv(last_time,rep(1,nrow(forsurv)))~new_alg,data=forsurv)
#For 5-modes with 100 simulations each
    # forsurv <- dist_t1_times  |> select(new_alg,last_time)
    # fit <- survfit(Surv(last_time,rep(1,nrow(forsurv)))~new_alg,data=forsurv)
#For 7-modes with 100 simulations each
    # forsurv <- dist_t1_times  |> select(new_alg,last_time)
    # fit <- survfit(Surv(last_time,rep(1,nrow(forsurv)))~new_alg,data=forsurv)
       
    forsurv <- dist_t1_times  |>
    select(new_alg,last_time) |> 
    mutate(status=ifelse(is.na(last_time),0,1))
    fit <- survfit(Surv(last_time,status)~new_alg,data=forsurv)
    
    (plot_surv_mode <- ggsurvplot(fit,
                                  data=forsurv,
                                  fun="event",
                                  palette = "Set1",    # Color palette
                                  xlab = "Time (seconds)",
                                  ylab = "Prop. of simulations visiting all modes",
                                  legend.title = "Algorithm",
                                  # break.time.by = time_br,
                                  font.x = 15,        # X-axis label font size
                                  font.y = 15,        # Y-axis label font size
                                  font.tickslab = 12, # Axis tick labels (numbers) font size
                                  font.legend = 10,
                                  conf.int = FALSE,
                                  censor = TRUE))   # Legend text font size)
    
    jpeg(file.path(export_path,paste0(export_file_name,"_speed_mode",".jpg")),width=1200,height =600,pointsize = 30)
    print(plot_surv_mode)
    dev.off()
    

    # ggplot(forsurv, aes(x = last_time, color = algorithm)) +
    #   stat_ecdf(size = 1) +  # Empirical CDF (inverse of survival)
    #   geom_point(
    #     data = filter(forsurv, status == 0),  # Censored points
    #     aes(y = 0.95), shape = 3, size = 3, stroke = 1.5
    #   ) +
    #   labs(
    #     y = "Prop. of simulations visiting all modes", 
    #     x = "Time (seconds)",
    #     #title = "Algorithm Performance (CDF)"
    #   ) +
    #   theme_minimal()
    
    
### Speed to modes for SECOND temperature    
    # dist_t2 <- distances |> 
    #   group_by(alg) |> 
    #   summarise(max_temp=max(temperature)) |> 
    #   ungroup() |> 
    #   right_join(distances,by="alg") |> 
    #   filter(temperature!=max_temp) |> 
    #   select(-sim,-max_temp) |> 
    #   group_by(alg) |> 
    #   summarise(max_temp=max(temperature)) |> 
    #   ungroup() |> 
    #   right_join(distances,by="alg") |> 
    #   filter(temperature==max_temp)
    # 
    # dist_t2_times <- dist_t2 |>  select(-min_dist) |> 
    #   pivot_longer(cols=c(time_find),names_to="variable",values_to="measure") |> 
    #   pivot_wider(names_from=mode,values_from=measure) |> 
    #   rowwise() |> 
    #   mutate(first_time=min(m1,m2),last_time=max(m1,m2))
    # 
    # 
    # forsurv <- dist_t2_times  |>
    #   select(algorithm,last_time)
    # 
    # fit <- survfit(Surv(last_time,rep(1,nrow(forsurv)))~algorithm,data=forsurv)
    # 
    # (plot_surv_mode <- ggsurvplot(fit,
    #                               data=forsurv,
    #                               fun="event",
    #                               palette = "Set1",    # Color palette
    #                               xlab = "Time (seconds)",
    #                               ylab = "Prop. of simulations visiting all modes",
    #                               legend.title = "Algorithm",
    #                               # break.time.by = time_br,
    #                               font.x = 15,        # X-axis label font size
    #                               font.y = 15,        # Y-axis label font size
    #                               font.tickslab = 12, # Axis tick labels (numbers) font size
    #                               font.legend = 10,
    #                               conf.int = FALSE))   # Legend text font size)
    


}




