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
# chosen_dim <- "highdim"; file_dim <- "highd"
chosen_dim <- "lowdim";file_dim <- "lowd" #,10000,1000,5000
print_bimodal <- FALSE
print_multimodal <- FALSE
chosen_ids <-c(182:185)
# chosen_ids <-c(130,132,134,136,138,142,143,144)#c(129,131,133,135,137,139,140,141)

#### Chosen for lowdim bimodal problem ####
#We stay with 3 temperatures for everything
#For PT-IIT the last temperature does not achieve the 0.23 swap rate (it gets bigger rate)
#For PT A-IITm there doesn't seem to be a big issue
#For PT A-IITw it seems to have been better with 4 temperatures
# but still doesn't get the needed swap rate

chosen_bimodal <- c(129,135,137,139)
print_bimodal <- TRUE
chosen_ids <- chosen_bimodal

#### Chosen for lowdim multimodal problem ####
# For PT-IIT and PT A-IITw we only use 3 temperatures because it had the best performance in TVD
# For PT A-IITm we still have to identify the best temperature
# 
# chosen_multimodal <- c(165,167,169, 190)
# print_multimodal <- TRUE
# chosen_ids <- chosen_multimodal


#List of files
parameters <- read_csv(paste0("results/simulation_details_",file_dim,".csv"))
#Create table with available files
data_sum <- tibble(file_names=list.files(path = "results", pattern = "^sim_.*\\.Rds")) |> 
  mutate(id=as.numeric(str_extract(file_names, "(?<=id_)[0-9]+(?=\\.Rds)")),
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
    round_trip <- as.data.frame(matrix(NA,ncol=(num_replicas+2))); colnames(round_trip) <- c("alg","sim",1:num_replicas)
    swap_rate <- as.data.frame(matrix(NA,ncol=(num_replicas+1))); colnames(swap_rate) <- c("alg","sim",1:(num_replicas-1))
    iterations <- as.data.frame(matrix(NA,ncol=(num_replicas+2))); colnames(iterations) <- c("alg","sim",1:num_replicas)
    
    list_of_states <- list()
    iter_visit <- as.data.frame(matrix(NA,ncol=max(data_sum$simulations)+2));colnames(iter_visit) <- c("alg","state",1:max(data_sum$simulations))
    loglik_visited <- as.data.frame(matrix(NA,ncol=4));colnames(loglik_visited) <- c("alg","sim","state","loglik")
    max_lik <- as.data.frame(matrix(NA,ncol=4));colnames(max_lik) <- c("alg","sim","iteration","lik")
  }

}
time_taken <- as.data.frame(matrix(NA,ncol=2));colnames(time_taken) <- c("alg","time")

full_iter <- list()
k <- 1
Q <- 1

# Start creating datasets with information
for(i in 1:nrow(data_sum)){
  data <- readRDS(file.path("results",data_sum[i,1]))
  tot_sim <- data_sum |> slice(i)|> pull(simulations)
  algorithm <- data_sum |> slice(i) |> pull(algorithm)
  if(algorithm=="PT_A_IIT"){algorithm <- "PT A-IIT m"}
  if(algorithm=="PT_IIT_no_Z"){algorithm <- "PT-IIT no Z"}
  if(algorithm=="PT_IIT_Z"){algorithm <- paste0("PT-IIT (",data_sum |> slice(i)|> pull(bf),")")}
  if(algorithm=="PT_A_IIT_RF"){algorithm <- "PT A-IIT w"}
  
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
  temp$sim <- 1:tot_sim
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
  temp$sim <- 1:tot_sim
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
    
    temp$sim <- 1:tot_sim
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
    temp$sim <- 1:tot_sim
    temp$alg <- algorithm
    temp <- temp |> select(alg,sim,everything())
    temp <- temp |> pivot_longer(cols=starts_with("V"),names_to="measurement",values_to="time")
    
    temp2 <- data[["tvd_report"]]
    temp2 <- as.data.frame(temp2)
    temp2$sim <- 1:tot_sim
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
    if(data_sum|> slice(i) |> pull(model)=="bimodal"){
      p <- data_sum$p;
    }
    if(data_sum|> slice(i) |> pull(model)=="gset"){
      file_matrix <- paste0("gset/",data_sum |> slice(i)|> pull(file),".txt")
      p <- readParameters(file_matrix)
      
 ### Evaluate maximum likelihood found     
      state_max <- matrix(NA,nrow=p,ncol=tot_sim)
      likelihood_max <- rep(0,tot_sim)
      #Extract specific state
      for(i in 1:tot_sim){#Extract 1 state for each simulation
        vec_temp <- data[["states"]][,index_max[i],i]
        state_max[,i] <- vec_temp
        # likelihood_max[i] <- eval_lik(file_matrix,vec_temp)
      }
      likelihood_max <- eval_lik_matrix(file_matrix,state_max)
      temp$lik <- likelihood_max
      max_lik <- rbind(max_lik,temp)
      
    }
    
# Each column is a simulation   
# Identifying maximum likelihood visited in each simulation
    index_max <- apply(data[["loglik_visited"]],2,which.max)
    full_index_max <- cbind(index_max,1:tot_sim)
    
    
    #Storing number of iterations to visit the state with max likelihood and the value of likelihood
    temp <- data[["iter_visit"]]
    temp <- as.data.frame(temp[full_index_max])
    colnames(temp) <- "iteration"
    temp$sim <- 1:tot_sim
    temp$alg <- algorithm
    temp <- temp |> select(alg,sim,everything())
    

    
    
# ###### Storing number of iterations to visit modes
#     temp <- as.data.frame(data[["iter_visit"]])
#     colnames(temp) <- 1:ncol(temp)
#     temp$state <- 1:nrow(temp)
#     temp$alg <- algorithm
#     temp <- temp |> select(alg,state,everything())
#     iter_visit <- rbind(iter_visit,temp)
###### Storing likelihood of visited states
    temp <- as.data.frame(data[["loglik_visited"]])
    #As many rows as number of states we were tracking
    number_sim <- ncol(temp)
    colnames(temp) <- paste0("s",1:number_sim) #As many columns as simulations
    temp$state <- 1:nrow(temp)
    temp$alg <- algorithm
    temp <- temp  |> select(alg,state,everything()) |> 
      pivot_longer(cols=all_of(paste0("s",1:number_sim)),names_to = "sim",values_to='loglik')
    
    #> 
    loglik_visited <- rbind(loglik_visited,temp)

    #Storing full list of states
    list_of_states[[Q]] <- data[["states"]]
    Q <- Q+1
  }
  
  
  if(!is_empty(data[["round_trips"]])){
    #Extract number of round trips rate
    temp <- as.data.frame(data[["round_trips"]])
    colnames(temp) <- 1:ncol(temp)
    temp$sim <- 1:tot_sim
    temp$alg <- algorithm
    temp <- temp |> select(alg,sim,everything())
    round_trip <- rbind(round_trip,temp)
    # Extract replica swap rate
    temp <- as.data.frame(data[["swap_rate"]])
    colnames(temp) <- 1:ncol(temp)
    temp$sim <- 1:tot_sim
    temp$alg <- algorithm
    temp <- temp |> select(alg,sim,everything())
    swap_rate <- rbind(swap_rate,temp)
  }
  if(!is_empty(data[["total_iter"]])){ 
    # Extract total iterations
    dims<- dim(data[["total_iter"]])
    full_iter[[k]] <- data[["total_iter"]]
    k <- k+1;
    temp <- as.data.frame(t(colSums(data[["total_iter"]])))/dims[1]
    #temp is the average number of Rejection Free steps before trying a swap
    colnames(temp) <- 1:ncol(temp)
    temp$sim <- 1:tot_sim
    temp$alg <- algorithm
    temp <- temp |> select(alg,sim,everything())
    iterations <- rbind(iterations,temp)
  }
}
##### Delete first row with NA#####
round_trip <- round_trip |> filter(!is.na(alg))
swap_rate <- swap_rate |> filter(!is.na(alg))
iterations <- iterations |> filter(!is.na(alg))
time_taken <- time_taken |> filter(!is.na(alg))


##### Export plots and tables #####
export_path <- paste0("C:/Users/ralex/Documents/src/paper_iit_pt/images/",chosen_dim,"_ex")
export_file_name <- paste0(chosen_dim,"_",paste0(chosen_ids,collapse="_"))
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

jpeg(file.path(export_path,paste0("time_table_",export_file_name,".jpg")),width=80*ncol(time_table),height=35*nrow(time_table),pointsize = 30)
grid.arrange(tableGrob(time_table))
dev.off()

(time_boxplot <- time_taken |> 
  filter(!str_starts(alg,'IIT')) |>
  ggplot(aes(x=alg,y=time))+
  geom_boxplot())

##### Report on average swap rate

#First create column names depending on the number of replicas
swap_names <- c()
for(i in 1:(num_replicas-1)){
  swap_names <- c(swap_names,paste0(i,"↔",i+1))
}

swap_report <- swap_rate |> 
  group_by(alg) |>
  summarize(across(-sim, function(x) mean(x, na.rm = TRUE)))
colnames(swap_report) <- c("alg",swap_names)

jpeg(file.path(export_path,paste0("table_swap_rate_",export_file_name,".jpg")),width=93*ncol(swap_report),height=33*nrow(swap_report),pointsize = 30)
grid.arrange(tableGrob(swap_report))
dev.off()

##### Report on average round trip rate
rt_report <- round_trip |> 
  group_by(alg) |>
  summarize(across(-sim, function(x) mean(x, na.rm = TRUE)))
colnames(rt_report) <- c("alg",paste0("R",1:num_replicas))

jpeg(file.path(export_path,paste0("table_roundtrip_rate_",export_file_name,".jpg")),width=90*ncol(rt_report),height=40*nrow(rt_report),pointsize = 30)
grid.arrange(tableGrob(rt_report))
dev.off()

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
  
  jpeg(file.path(export_path,paste0("tvd_",export_file_name,".jpg")),width=800,height =400,pointsize = 30)
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

jpeg(file.path(export_path,paste0("table_visited_modes_",export_file_name,".jpg")),width=100 + (60*(ncol(table_visited)-1)),height=40*nrow(table_visited),pointsize = 30)
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

jpeg(file.path(export_path,paste0("table_iterations_",export_file_name,".jpg")),width=150*nrow(iterations_to_explore),height=40*nrow(iterations_to_explore),pointsize = 30)
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
  
}

jpeg(file.path(export_path,paste0("table_swaps_",export_file_name,".jpg")),width=140*nrow(swaps_to_explore),height=40*nrow(swaps_to_explore),pointsize = 30)
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

fit <- survfit(Surv(last_time,rep(1,nrow(forsurv)))~alg,data=forsurv)
if(print_bimodal){time_br <- 0.2}
if(print_multimodal){time_br <- 0.5}
(plot_surv_mode <- ggsurvplot(fit,
           data=forsurv,
           fun="event",
           palette = "Set1",    # Color palette
           xlab = "Time (seconds)",
           ylab = "Prop. of simulations visiting all modes",
           legend.title = "Algorithm",
           break.time.by = time_br,
           font.x = 15,        # X-axis label font size
           font.y = 15,        # Y-axis label font size
           font.tickslab = 12, # Axis tick labels (numbers) font size
           font.legend = 16))   # Legend text font size)

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

jpeg(file.path(export_path,paste0("time_mode_",export_file_name,".jpg")),width=800,height =400,pointsize = 30)
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



jpeg(file.path(export_path,paste0("tvd_mean_",export_file_name,".jpg")),width=800,height =400,pointsize = 30)
print(mean_tvd_plot)
dev.off()

jpeg(file.path(export_path,paste0("tvd_max_",export_file_name,".jpg")),width=800,height =400,pointsize = 30)
print(max_tvd_plot)
dev.off()


}# Finish low dim reports

#Starts high dim reports
if(chosen_dim=="highdim"){
  ##### Delete first row with NA#####
  max_lik <- max_lik|> filter(!is.na(alg))
  loglik_visited <- loglik_visited |> filter(!is.na(alg))
  iter_visit <- iter_visit|> filter(!is.na(alg))

}




