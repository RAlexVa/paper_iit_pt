rm(list=ls())
### Process simulation results
library(stringr)
library(tidyverse)
# data_sum <- data.frame(file=dir("results"))
# data_sum$alg <- str_extract(results, ".*(?=_sim)")


data_sum <- tibble(file = dir("results")) |> 
  mutate(alg = str_extract(file, ".*(?=_sim)"),
         sim = str_extract(file,"(?<=sim_)[^_]+"),
         iter=str_extract(file,"(?<=iter_)[^(_\\.)]+"),
         interswap=str_extract(file,"(?<=interswap_)[^_]+"),
         tot_swap=str_extract(file,"(?<=totalswap_)[^\\.]+"),
         iterswap=str_extract(file,"(?<=iterswap_)[^\\.]+")) |> 
  mutate(across(-(1:2), as.numeric))

data_sum <- data_sum |> 
  mutate(interswap=ifelse(is.na(interswap),iterswap,interswap)) |> 
  mutate(iter=if_else(is.na(iter),interswap*(tot_swap),iter)) |> 
  mutate(tot_swap=if_else(is.na(tot_swap),iter/interswap,tot_swap)) |> 
  select(-iterswap)

tvd <- data.frame(alg=NA,sim=NA,tvd=NA)
mode_visit <- as.data.frame(matrix(NA,ncol=9)); colnames(mode_visit) <- c("alg","sim",1:7)
round_trip <- as.data.frame(matrix(NA,ncol=6)); colnames(round_trip) <- c("alg","sim",1:4)
swap_rate <- as.data.frame(matrix(NA,ncol=5)); colnames(swap_rate) <- c("alg","sim",1:3)
iterations <- as.data.frame(matrix(NA,ncol=6)); colnames(iterations) <- c("alg","sim",1:4)
for(i in 1:nrow(data_sum)){
  data <- readRDS(file.path("results",data_sum[i,1]))
  tot_sim <- data_sum |> slice(i) |> pull(sim)
  algorithm <- data_sum |> slice(i) |> pull(alg)
  print(data_sum[i,"alg"])
  print(names(data))

# Extract TVD
  temp <- tibble(alg=algorithm,sim=1:tot_sim,tvd=data[["tvd"]])
  tvd <- rbind(tvd,temp)
  
# Extract visit of modes
  temp <- as.data.frame(data[["mode_visit"]])
  colnames(temp) <- 1:ncol(temp)
  temp$sim <- 1:tot_sim
  temp$alg <- algorithm
  temp <- temp |> select(alg,sim,everything())
  mode_visit <- rbind(mode_visit,temp)

  if(algorithm!='IIT'){
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
  if(algorithm=="PT_A_IIT"){ 
    # Extract total iterations
    dims<- dim(data[["total_iter"]])
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
tvd <-  tvd |> filter(!is.na(alg))
mode_visit <- mode_visit |> filter(!is.na(alg))
round_trip <- round_trip |> filter(!is.na(alg))
swap_rate <- swap_rate |> filter(!is.na(alg))
iterations <- iterations |> filter(!is.na(alg))
##### Total Variation Distance #####

tvd |>  filter(alg!='IIT') |> 
  ggplot(aes(x=alg,y=tvd,fill=alg)) +
  geom_boxplot()+
  labs(fill='Algortihm',x="",y="Total Variation Distance")

##### First visit to modes #####
col_names <- c("1","2","3","4","5","6","7")

mode_sum <- mode_visit |> 
  rowwise()|> 
  mutate(last_visit=max(c_across(col_names)), first_visit=min(c_across(col_names)[c_across(col_names)>0])) |> 
  mutate(first_mode =  names(pick(col_names))[which(c_across(col_names) == first_visit)[1]]) |> 
  mutate(last_mode =  names(pick(col_names))[which(c_across(col_names) == last_visit)[1]]) |> 
  mutate(total_modes = sum(c_across(col_names) > 0))
  



