rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(dplyr)
library(readr)
# install.packages("Rcpp")
# install.packages("RcppArmadillo")
# install.packages("readr")
# install.packages("tidyverse")

##### import functions #####

source("functions/r_functions.R")
run_highd <- function(list_ids){
  if(!("./results" %in% list.dirs(recursive=F))){
    print("Wrong directory. There's no results folder for the output")
  }else{
    
  if(missing(list_ids)){
    #### Prompt to choose which simulation to run
    writeLines("You can write various IDs separated by commas")
    list_ids <- readline('Choose id:')
  }
    if(!is.character(list_ids)){ list_ids <- as.character(list_ids)}
    list_ids <- as.numeric(unlist(strsplit(list_ids,",")))
  ##### Read file for parameters #####
  parameters <- as.data.frame(read_csv("results/simulation_details_highd.csv"))
  
#Check how many models we're doing
    tot_models <- unique(parameters|> filter(id %in% list_ids) |> pull(model))
    if(length(tot_models)==1){only_1_model <- TRUE;}else{only_1_model <- FALSE}
# In case only 1 model is chosen, we only read 1 model for all the IDs    
    if(only_1_model){
      print("Reading one set of C++ functions")
      if(tot_models=="gset"){Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp")}
      if(tot_models=="bimodal"){Rcpp::sourceCpp("functions/cpp_func_multihigh.cpp")}
    }
  #Start process for algorithms
  for(id_chosen in list_ids){

    sim_chosen <- parameters |> filter(id==id_chosen)
    if(nrow(sim_chosen)!=1){print(paste0("Error: id ",id_chosen," doesn't exist or there's more than one")); next;}
#In case there are more than 1 model, I need to re-read functions depending on model
    if(!only_1_model){
      print(paste0("Reading C++ functions for id: ",id_chosen," model: ",sim_chosen$model))
      if(sim_chosen$model=="gset"){Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp")}
      if(sim_chosen$model=="bimodal"){Rcpp::sourceCpp("functions/cpp_func_multihigh.cpp")}
    }
    # Parameters for all algorithms
    total_simulations <- sim_chosen$simulations
    temperatures <- as.numeric(sim_chosen |> select(matches("^t\\d+$")))
    bal_f <- as.character(sim_chosen|> select(matches("^bf")))
    defined_seed <- sim_chosen$seed
    set.seed(defined_seed)
    #Parameters for PT with IIT
    total_iter <- sim_chosen$iterations #300000 #Total number of steps to perform in each replica
    iterswap <- sim_chosen$interswap #Total iterations before trying a replica swap
    #Parameters for PT with a-IIT
    sample_inter_swap <- sim_chosen$interswap #Number of original samples to get before trying a replica swap
    total_swap <- sim_chosen$total_swap #Total number of swaps to try
    burnin_iter <- sim_chosen$burn_in #Number of iterations for burn-in
    file_matrix <- paste0("gset/",sim_chosen$file,".txt")
    p <- readParameters(file_matrix)
    states_visited <- sim_chosen$states_visited
    # start_state <- sim_chosen$start_state;
    alg <- sim_chosen$algorithm
    start_state <- as.numeric(unlist(strsplit(sim_chosen$start_state,",")))
    export <- list();
    #### Function depending on algorithm to use
    
    writeLines(c("Parameters:",paste0("Algorithm: ",alg),
                 paste0("ID: ",id_chosen),
                 paste0("Seed: ",defined_seed),
                 paste0("Total simulations: ",total_simulations),
                 paste0("Burn-in iterations: ",burnin_iter),
                 paste0("Temperatures: ",paste(temperatures,collapse=',')),
                 paste0("Balancing functions: ",paste(bal_f,collapse = ',')),
                 paste0("Total iterations: ",total_iter),
                 paste0("Try swaps:",iterswap),
                 paste0("Samples in-between swaps: ",sample_inter_swap),
                 paste0("Total swaps:",total_swap),
                 paste0("File: ",file_matrix),
                 paste0("States to keep track: ",states_visited)))
    
    # check <- as.numeric(readline('ok? 1 Yes/ 0 No'))
    check <- 1;
    if(check!=1){print("modify parameters")}else{
      if(alg=="IIT"){
        # Only IIT
        # PT_IIT_sim(int p,int startsim,int endsim, int numiter, int iterswap,int burn_in, vec temp, const std::vector<std::string>& bal_function, bool bias_fix,const std::string& filename,int num_states_visited)
        output <- PT_IIT_sim(p,1,total_simulations,total_iter,total_iter+1,burnin_iter,temperatures[1],bal_f[1], TRUE, file_matrix,states_visited,start_state)
      }else{
        
        if(alg=="PT_IIT_Z"){
          # Using Z factor bias correction
          #PT_IIT_sim(int p,int startsim,int endsim, int numiter, int iterswap,int burn_in, vec temp, const std::vector<std::string>& bal_function, bool bias_fix,const std::string& filename,int num_states_visited)
          output <- PT_IIT_sim(p,1,total_simulations,total_iter,iterswap,burnin_iter,temperatures,bal_f,TRUE, file_matrix,states_visited,start_state)
          #round trip rate (NA for IIT)
          export[["round_trips"]] <- PT_RT(output[["ip"]], floor(total_iter/iterswap),total_simulations)
        }
        if(alg=="PT_IIT_no_Z"){
          # output_name <- paste0("PT_IIT_no_Z_","sim_",total_simulations,"_iter_",total_iter,"_iterswap_",iterswap,"_s_",defined_seed,".Rds");
          # Without Z factor bias correction
          #PT_IIT_sim(int p,int startsim,int endsim, int numiter, int iterswap,int burn_in, vec temp, const std::vector<std::string>& bal_function, bool bias_fix,const std::string& filename,int num_states_visited)
          output <- PT_IIT_sim(p,1, total_simulations,total_iter,iterswap,burnin_iter,temperatures,bal_f,FALSE, file_matrix,states_visited,start_state)
          #round trip rate (NA for IIT)
          export[["round_trips"]] <- PT_RT(output[["ip"]], floor(total_iter/iterswap),total_simulations)
        }
        if(alg=="PT_A_IIT"){
          # Using A-IIT in each replica
          #PT_a_IIT_sim(int p,int startsim,int endsim, int total_swaps,int sample_inter_swap,int burn_in, vec temp, const std::vector<std::string>& bal_function,const std::string& filename,int num_states_visited)
          output <- PT_a_IIT_sim(p,1,total_simulations,total_swap,sample_inter_swap,burnin_iter,temperatures,bal_f,file_matrix,states_visited,start_state)
          #Number of iterations needed between swaps for each replica
          export[["total_iter"]] <- output[["total_iter"]]
          #round trip rate (NA for IIT)
          export[["round_trips"]] <- PT_RT(output[["ip"]],total_swap,total_simulations)
        }
        if(alg=="PT_A_IIT_RF"){
          # Using A-IIT with weights in each replica
          #PT_a_IIT_sim_RF(int p,int startsim,int endsim, int numiter, int iterswap,int burn_in, vec temp, const std::vector<std::string>& bal_function, bool bias_fix,const std::string& filename,int num_states_visited,const std::vector<int>& starting_coord)
          output <- PT_a_IIT_sim_RF(p,1,total_simulations,total_iter,iterswap,burnin_iter,temperatures,bal_f,TRUE,file_matrix,states_visited,start_state)
          #round trip rate (NA for IIT)
          export[["round_trips"]] <- PT_RT(output[["ip"]], floor(total_iter/iterswap),total_simulations)
        }
        # Replica swap acceptance rate (NA for IIT)
        export[["swap_rate"]] <- output[["swap_rate"]]
      }
        # Replica swap acceptance rate (NA for IIT)
        export[["swap_rate"]] <- output[["swap_rate"]]
        
      }
      
      export[["states"]] <- output[["states"]]
      export[["loglik_visited"]] <- output[["loglik_visited"]]
      export[["iter_visit"]]<- output[["iter_visit"]]
      export[["time_taken"]] <- output[["time_taken"]]
      export[["ip"]] <- output[["ip"]]
      output_name <- paste0("sim_highdim_id_",id_chosen,".Rds")
      saveRDS(export,file=file.path("results",output_name))
    }
  }
}

