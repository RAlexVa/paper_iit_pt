rm(list=ls())
# Import libraries
library(stringr)
library(tidyverse)
library(gridExtra)# For tables
library(latex2exp) #For using latex

### Process simulation results ###
#Define seeds to filter
filtered_seed <- 153
#Create table with available files
data_sum <- tibble(file = dir("results")) |> 
  mutate(alg = str_extract(file, ".*(?=_sim)"),
         sim = str_extract(file,"(?<=sim_)[^_]+"),
         iter=str_extract(file,"(?<=iter_)[^_]+"),
         interswap=str_extract(file,"(?<=interswap_)[^_]+"),
         tot_swap=str_extract(file,"(?<=totalswap_)[^_]+"),
         iterswap=str_extract(file,"(?<=iterswap_)[^_]+"),
         seed=str_extract(file,"(?<=s_)[^\\.]+")) |> 
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
full_iter <- list()
k <- 1

# filter specified seed
 data_sum <- data_sum |> filter(seed==filtered_seed)
# Start creating datasets with information
for(i in 1:nrow(data_sum)){
  data <- readRDS(file.path("results",data_sum[i,1]))
  tot_sim <- data_sum |> slice(i) |> pull(sim)
  algorithm <- data_sum |> slice(i) |> pull(alg)
  if(algorithm=="PT_A_IIT"){algorithm <- "A-IIT"}
  if(algorithm=="PT_IIT_no_Z"){algorithm <- "IIT no Z"}
  if(algorithm=="PT_IIT_Z"){algorithm <- "IIT w Z"}
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
  if(algorithm=="A-IIT"){ 
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
tvd <-  tvd |> filter(!is.na(alg))
mode_visit <- mode_visit |> filter(!is.na(alg))
round_trip <- round_trip |> filter(!is.na(alg))
swap_rate <- swap_rate |> filter(!is.na(alg))
iterations <- iterations |> filter(!is.na(alg))

##### Export plots and tables #####
export_path <- "C:/Users/ralex/Documents/src/paper-adaptive-iit-latex/images/lowdim_ex"
##### Total Variation Distance #####

tvd_plot <- tvd |>  filter(alg!='IIT') |> 
  ggplot(aes(x=alg,y=tvd,fill=alg)) +
  geom_boxplot(show.legend = FALSE)+
  labs(fill='Algortihm',x="",y="Total Variation Distance")+
  theme_minimal(base_size = 17)+
  theme(legend.key.size = unit(1, 'cm'))
tvd_plot

jpeg(file.path(export_path,"tvd_lowdim.jpg"),width=800,height =400,pointsize = 30)
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

jpeg(file.path(export_path,"table_visited_modes.jpg"),,width=50*nrow(table_visited),height=40*nrow(table_visited),pointsize = 30)
grid.arrange(tableGrob(table_visited))
dev.off()




##### Report on number of iterations for the original replica to visit all modes (most of the times after a swap)
iterations_to_explore <- mode_sum |> filter(alg!='IIT') |> 
  group_by(alg) |> 
  summarise(min=min(last_visit),
            q1=quantile(last_visit,probs=0.25),
            median=quantile(last_visit,probs=0.5),
            mean=mean(last_visit),
            q3=quantile(last_visit,probs=0.75),
            max=max(last_visit))

jpeg(file.path(export_path,"table_iterations.jpg"),width=150*nrow(iterations_to_explore),height=40*nrow(iterations_to_explore),pointsize = 30)
grid.arrange(tableGrob(iterations_to_explore))
dev.off()


##### Report on number of replica swaps needed to visit all modes
# First algorithms do 2k iterations before trying a replica swap
swaps_to_explore <- mode_sum |> filter(alg!='IIT',alg!="A-IIT") |> 
  group_by(alg) |> 
  summarise(min=min(last_visit)/2000,
            q1=quantile(last_visit,probs=0.25)/2000,
            median=quantile(last_visit,probs=0.5)/2000,
            mean=mean(last_visit)/2000,
            q3=quantile(last_visit,probs=0.75)/2000,
            max=max(last_visit)/2000)

temp <- mode_sum |> filter(alg=="A-IIT") |>
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
jpeg(file.path(export_path,"table_swaps.jpg"),,width=140*nrow(swaps_to_explore),height=40*nrow(swaps_to_explore),pointsize = 30)
grid.arrange(tableGrob(swaps_to_explore))
dev.off()


##### Report on average swap rate
swap_report <- swap_rate |> 
  group_by(alg) |>
  summarise(`1↔2`=mean(`1`),`2↔3`=mean(`2`),`3↔4`=mean(`3`))


jpeg(file.path(export_path,"table_swap_rate.jpg"),,width=140*nrow(swap_report),height=40*nrow(swap_report),pointsize = 30)
grid.arrange(tableGrob(swap_report))
dev.off()


##### Report on average swap rate
rt_report <- round_trip |> 
  group_by(alg) |>
  summarise(`R1`=mean(`1`),`R2`=mean(`2`),`R3`=mean(`3`),,`R4`=mean(`4`))

jpeg(file.path(export_path,"table_roundtrip_rate.jpg"),,width=140*nrow(rt_report),height=40*nrow(rt_report),pointsize = 30)
grid.arrange(tableGrob(rt_report))
dev.off()


