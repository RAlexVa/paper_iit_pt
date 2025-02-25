rm(list=ls())
# Import libraries
library(stringr)
library(tidyverse)
library(gridExtra)# For tables
library(latex2exp) #For using latex


### Process simulation results ###

# Choose dimension
# chosen_dim <- "highdim"; file_dim <- "highd"
chosen_dim <- "lowdim";file_dim <- "lowd" #,10000,1000,5000
chosen_ids <-c(21,22,23,24)#c(28,29,30,31,32,33)#c(31,32)#c(29,30)#c(25,26,27)#c(17,18,19,20)#c(25,26,27)#c(9,10,11,12)   #c(20,21,22) # c(13,14,15,16)



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
  #Low dim is the example with 7 modes and 4 temperatures
  tvd <- data.frame(alg=NA,sim=NA,tvd=NA)
  mode_visit <- as.data.frame(matrix(NA,ncol=10)); colnames(mode_visit) <- c("alg","sim","interswap",1:7)
  round_trip <- as.data.frame(matrix(NA,ncol=6)); colnames(round_trip) <- c("alg","sim",1:4)
  swap_rate <- as.data.frame(matrix(NA,ncol=5)); colnames(swap_rate) <- c("alg","sim",1:3)
  iterations <- as.data.frame(matrix(NA,ncol=6)); colnames(iterations) <- c("alg","sim",1:4) #For the 4 replicas
  pi_modes <- as.data.frame(matrix(NA,ncol=9));colnames(pi_modes) <- c("alg","sim",1:7)
  # Low dimensional true probability setup
  {
    Rcpp::sourceCpp("functions/cpp_functions.cpp") #To use vec_to_num function
    source("functions/r_functions.R")
    p <- 16 #dimension
    theta <- 8 #tail weight parameter
    
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
      pi.true[i+1] <-  ll_comp(NumberToVector(i,p),modes_list,theta,p)
    }
    pi.true <- exp(pi.true)/(length(modes_list)*(1+exp(-theta))^p)
  }

  
  
}
if(chosen_dim=="highdim"){
  Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp") #To use eval_lik function
  source("functions/r_functions.R") #To get dimension of the problem
  #High dim we have 3 temperatures
  #we dont track distribution convergence
  #we track visit to high probability states
  round_trip <- as.data.frame(matrix(NA,ncol=5)); colnames(round_trip) <- c("alg","sim",1:3)
  swap_rate <- as.data.frame(matrix(NA,ncol=4)); colnames(swap_rate) <- c("alg","sim",1:2)
  iterations <- as.data.frame(matrix(NA,ncol=5)); colnames(iterations) <- c("alg","sim",1:3)
  
  list_of_states <- list()
  iter_visit <- as.data.frame(matrix(NA,ncol=max(data_sum$simulations)+2));colnames(iter_visit) <- c("alg","state",1:max(data_sum$simulations))
  loglik_visited <- as.data.frame(matrix(NA,ncol=max(data_sum$simulations)+2));colnames(loglik_visited) <- c("alg","state",1:max(data_sum$simulations))
  max_lik <- as.data.frame(matrix(NA,ncol=4));colnames(max_lik) <- c("alg","sim","iteration","lik")
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
  if(algorithm=="PT_IIT_Z"){algorithm <- "PT-IIT"}
  if(algorithm=="PT_A_IIT_RF"){algorithm <- "PT A-IIT w"}
  
##### Optinal add the ID of the simulation into the name of the algorithm
  algorithm <- paste0(algorithm,"(",data_sum |> slice(i) |> pull(id),")")
  
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

}
  if(chosen_dim=="highdim"){
### Specific extractions for highdim example
    file_matrix <- paste0("gset/",data_sum |> slice(i)|> pull(file),".txt")
    p <- readParameters(file_matrix)
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
    
    
# ###### Storing number of iterations to visit modes
#     temp <- as.data.frame(data[["iter_visit"]])
#     colnames(temp) <- 1:ncol(temp)
#     temp$state <- 1:nrow(temp)
#     temp$alg <- algorithm
#     temp <- temp |> select(alg,state,everything())
#     iter_visit <- rbind(iter_visit,temp)
# ###### Storing likelihood of visited states
#     temp <- as.data.frame(data[["loglik_visited"]])
#     colnames(temp) <- 1:ncol(temp)
#     temp$state <- 1:nrow(temp)
#     temp$alg <- algorithm
#     temp <- temp |> select(alg,state,everything())
#     loglik_visited <- rbind(loglik_visited,temp)

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
# full_path <- file.path(export_path,paste0("tvd_",export_file_name,".jpg"))


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

jpeg(file.path(export_path,paste0("time_table_",export_file_name,".jpg")),width=250*nrow(time_table),height=35*nrow(time_table),pointsize = 30)
grid.arrange(tableGrob(time_table))
dev.off()



if(chosen_dim=="lowdim"){
  ##### Delete first row with NA#####
  tvd <-  tvd |> filter(!is.na(alg)) 
  mode_visit <- mode_visit |> filter(!is.na(alg))
  pi_modes <- pi_modes|> filter(!is.na(alg))
##### Compare estimation of modes #####  
  # pi_modes |> pivot_longer(cols = -(alg:sim), names_to = "mode", values_to = "pi_est") |> 
  #   ggplot(aes(x=alg,y=pi_est,fill=alg))+
  #   geom_boxplot(show.legend = FALSE)+
  #   facet_wrap(~mode)
  
  pi_modes |> pivot_longer(cols = -(alg:sim), names_to = "mode", values_to = "pi_est") |> 
    ggplot(aes(x=mode,y=pi_est,fill=alg))+
    geom_boxplot(show.legend = FALSE)+
    geom_hline(yintercept = pi.true[modes[1]+1], color = "red", linetype = "dashed", size = 1)+
    facet_wrap(~alg)+
    theme_minimal(base_size = 17)+
    theme(legend.key.size = unit(1, 'cm'))
  
  
  pi_modes |> rowwise() |> 
    mutate(tvd=0.5*sum(abs(pi.true[modes[1]+1]-c_across(`1`:`7`)))) |> 
    ungroup() |> 
    select(alg,tvd) |> 
    ggplot(aes(x=alg,y=tvd,fill=alg)) +
    geom_boxplot(show.legend = FALSE)+
    labs(fill='Algortihm',x="",y="Total Variation Distance")+
    theme_minimal(base_size = 17)+
    theme(legend.key.size = unit(1, 'cm'))
  
##### Total Variation Distance #####
  
  tvd_plot <- tvd |>  filter(!str_starts(alg,'IIT')) |> 
    ggplot(aes(x=alg,y=tvd,fill=alg)) +
    geom_boxplot(show.legend = FALSE)+
    labs(fill='Algortihm',x="",y="Total Variation Distance")+
    theme_minimal(base_size = 17)+
    theme(legend.key.size = unit(1, 'cm'))
  tvd_plot
  
  jpeg(file.path(export_path,paste0("tvd_",export_file_name,".jpg")),width=800,height =400,pointsize = 30)
  tvd_plot
  dev.off()
##### First visit to modes #####
col_names <- c("1","2","3","4","5","6","7")

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

jpeg(file.path(export_path,paste0("table_visited_modes_",export_file_name,".jpg")),width=50*nrow(table_visited),height=40*nrow(table_visited),pointsize = 30)
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

jpeg(file.path(export_path,paste0("table_swaps_",export_file_name,".jpg")),width=140*nrow(swaps_to_explore),height=40*nrow(swaps_to_explore),pointsize = 30)
grid.arrange(tableGrob(swaps_to_explore))
dev.off()

##### Report on average swap rate
swap_report <- swap_rate |> 
  group_by(alg) |>
  summarise(`1↔2`=mean(`1`),`2↔3`=mean(`2`),`3↔4`=mean(`3`))


jpeg(file.path(export_path,paste0("table_swap_rate_",export_file_name,".jpg")),width=140*nrow(swap_report),height=40*nrow(swap_report),pointsize = 30)
grid.arrange(tableGrob(swap_report))
dev.off()
##### Report on average swap rate
rt_report <- round_trip |> 
  group_by(alg) |>
  summarise(`R1`=mean(`1`),`R2`=mean(`2`),`R3`=mean(`3`),,`R4`=mean(`4`))

jpeg(file.path(export_path,paste0("table_roundtrip_rate_",export_file_name,".jpg")),width=140*nrow(rt_report),height=40*nrow(rt_report),pointsize = 30)
grid.arrange(tableGrob(rt_report))
dev.off()
}# Finish low dim reports

#Starts high dim reports
if(chosen_dim=="highdim"){
  ##### Delete first row with NA#####
  max_lik <- max_lik|> filter(!is.na(alg))
  loglik_visited <- loglik_visited |> filter(!is.na(alg))
  iter_visit <- iter_visit|> filter(!is.na(alg))
  ##### Report on average swap rate
  swap_report <- swap_rate |> 
    group_by(alg) |>
    summarise(`1↔2`=mean(`1`),`2↔3`=mean(`2`))
  
  
  jpeg(file.path(export_path,paste0("table_swap_rate_",export_file_name,".jpg")),width=140*nrow(swap_report),height=40*nrow(swap_report),pointsize = 30)
  grid.arrange(tableGrob(swap_report))
  dev.off()
  ##### Report on average swap rate
  rt_report <- round_trip |> 
    group_by(alg) |>
    summarise(`R1`=mean(`1`),`R2`=mean(`2`),`R3`=mean(`3`),)
  
  jpeg(file.path(export_path,paste0("table_roundtrip_rate_",export_file_name,".jpg")),width=140*nrow(rt_report),height=40*nrow(rt_report),pointsize = 30)
  grid.arrange(tableGrob(rt_report))
  dev.off()
}




