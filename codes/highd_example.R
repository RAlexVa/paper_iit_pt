rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(dplyr, warn.conflicts = F)
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
  parameters <- as.data.frame(read_csv("results/simulation_details_highd.csv", col_types = cols()))
  
#Check how many models we're doing
    tot_models <- unique(parameters|> filter(id %in% list_ids) |> pull(model))
    if(length(tot_models)==1){only_1_model <- TRUE;}else{only_1_model <- FALSE}
# In case only 1 model is chosen, we only read 1 model for all the IDs    
    if(only_1_model){
      print("Reading one set of C++ functions")
      print(paste0("Model = ",tot_models))
      if(tot_models=="gset"){Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp");
      file_matrix <- paste0("gset/",unique(parameters|> filter(id %in% list_ids) |> pull(file)),".txt");
      p <- readParameters(file_matrix);}
      if(tot_models=="bimodal"){Rcpp::sourceCpp("functions/cpp_functions_highdim_2.cpp");
        p <- unique(parameters|> filter(id %in% list_ids) |> pull(p));
        file_matrix <- "NA";
        if(length(p)>1){stop("You're trying to run multiple values for p")}}
    }
  #Start process for algorithms
  for(id_chosen in list_ids){

    sim_chosen <- parameters |> filter(id==id_chosen)
    if(nrow(sim_chosen)!=1){print(paste0("Error: id ",id_chosen," doesn't exist or there's more than one")); next;}
#In case there are more than 1 model, I need to re-read functions depending on model
    if(!only_1_model){
      print(paste0("Reading C++ functions for id: ",id_chosen," model: ",sim_chosen$model))
      if(sim_chosen$model=="gset"){Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp");
        file_matrix <- paste0("gset/",sim_chosen$file,".txt");
      p <- readParameters(file_matrix);}
      if(sim_chosen$model=="bimodal"){Rcpp::sourceCpp("functions/cpp_functions_highdim_2.cpp");
        p <- unique(parameters|> filter(id==id_chosen) |> pull(p));file_matrix <- "NA";}
    }
    # Parameters for all algorithms
    total_simulations <- sim_chosen$simulations
    temperatures <- as.numeric(sim_chosen |> select(matches("^t\\d+$")))
    temperatures <- temperatures[!is.na(temperatures)] #Ignore NA temperatures
    bal_f <- as.character(sim_chosen|> select(matches("^bf")))
    bal_f <- bal_f[!is.na(bal_f)] #Ignore NA balancing functions
    if(length(bal_f)==0){bal_f <- "sq"}
    burnin_iter <- sim_chosen$burn_in #Number of iterations for burn-in
    states_visited <- sim_chosen$states_visited# start state is random
    start_state <- -1
    alg <- sim_chosen$algorithm
    defined_seed <- sim_chosen$seed
    set.seed(defined_seed)
    #Parameters for PT with IIT
    total_iter <- sim_chosen$iterations 
    iterswap <- sim_chosen$interswap #Total iterations before trying a replica swap
    #Parameters for PT with A-IIT
    sample_inter_swap <- sim_chosen$interswap #Number of original samples to get before trying a replica swap
    total_swap <- sim_chosen$total_swap #Total number of swaps to try
    reduc_constant <- sim_chosen$reduc_constant;if(is.null(reduc_constant)){reduc_constant <- 0}
    reduc_model <- sim_chosen$reduc_model;if(is.null(reduc_model)){reduc_model <- "never"}
    
    export <- list();
    #### Function depending on algorithm to use
    
    writeLines(c("Parameters:",
                 paste0("ID: ",id_chosen),
                 paste0("Algorithm: ",alg),
                 paste0("Problem: ",sim_chosen$model),
                 paste0("Number of Replicas: ",length(temperatures)),
                 paste0("Total simulations: ",total_simulations),
                 paste0("Burn-in iterations: ",burnin_iter),
                 paste0("Total iterations: ",total_iter),
                 paste0("Temperatures: ",paste(temperatures,collapse=',')),
                 paste0("Balancing function: ",paste(bal_f,collapse = ',')),
                 paste0("Try swaps:",iterswap),
                 paste0("Samples in-between swaps: ",sample_inter_swap),
                 paste0("Total swaps:",total_swap),
                 paste0("File: ",file_matrix),
                 paste0("States to keep track: ",states_visited),
                 paste0("Reduction constant:",reduc_constant),
                 paste0("bound reduction method:",reduc_model)))
    
    #re-create the balancing function vector
    #All replicas use the same balancing function
    bal_f <- rep(bal_f,length(temperatures))
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
          #PT_IIT_sim(int p,int startsim,int endsim, int numiter, int iterswap,int burn_in, vec temp, const std::vector<std::string>& bal_function, bool bias_fix,const std::string& filename,int num_states_visited,const std::vector<int>& starting_coord)
          output <- PT_IIT_sim(p,1,total_simulations,total_iter,iterswap,burnin_iter,temperatures,bal_f,TRUE, file_matrix,states_visited,start_state)
          #round trip rate (NA for IIT)
          swaps_for_rt_rate <- floor(total_iter/iterswap)
        }
        if(alg=="PT_IIT_no_Z"){
          # output_name <- paste0("PT_IIT_no_Z_","sim_",total_simulations,"_iter_",total_iter,"_iterswap_",iterswap,"_s_",defined_seed,".Rds");
          # Without Z factor bias correction
          #PT_IIT_sim(int p,int startsim,int endsim, int numiter, int iterswap,int burn_in, vec temp, const std::vector<std::string>& bal_function, bool bias_fix,const std::string& filename,int num_states_visited,const std::vector<int>& starting_coord)
          output <- PT_IIT_sim(p,1, total_simulations,total_iter,iterswap,burnin_iter,temperatures,bal_f,FALSE, file_matrix,states_visited,start_state)
          #round trip rate (NA for IIT)
          swaps_for_rt_rate <- floor(total_iter/iterswap)
        }
        if(alg=="PT_A_IIT"){
          # Using A-IIT in each replica
          #PT_a_IIT_sim(int p,int startsim,int endsim, int total_swaps,int sample_inter_swap,int burn_in, vec temp, const std::vector<std::string>& bal_function,const std::string& filename,int num_states_visited,const std::vector<int>& starting_coord, double decreasing_constant,std::string reduc_model)
          output <- PT_a_IIT_sim(p,1,total_simulations,total_swap,sample_inter_swap,burnin_iter,temperatures,bal_f,file_matrix,states_visited,start_state,reduc_constant,reduc_model)
          #round trip rate (NA for IIT)
          swaps_for_rt_rate <- total_swap
        }
        if(alg=="PT_A_IIT_RF"){
          # Using A-IIT with weights in each replica
          #PT_a_IIT_sim_RF(int p,int startsim,int endsim, int numiter, int iterswap,int burn_in, vec temp, const std::vector<std::string>& bal_function, bool bias_fix,const std::string& filename,int num_states_visited,const std::vector<int>& starting_coord, double decreasing_constant,std::string reduc_model)
          output <- PT_a_IIT_sim_RF(p,1,total_simulations,total_iter,iterswap,burnin_iter,temperatures,bal_f,TRUE,file_matrix,states_visited,start_state,reduc_constant,reduc_model)
          #round trip rate (NA for IIT)
          swaps_for_rt_rate <- floor(total_iter/iterswap)
        }
      }
      #First export everything
      saveRDS(output,file=file.path("results",paste0("raw_sim_highdim_id_",id_chosen,".Rds")))
      for(e in names(output)){
        export[[e]] <- output[[e]]
      }
      #Then the exports that depend on the algorithm
      if(!is.null(output[["ip"]])){
        export[["round_trips"]] <- PT_RT(output[["ip"]], swaps_for_rt_rate,total_simulations)
      }
        output_name <- paste0("sim_highdim_id_",id_chosen,".Rds")
        saveRDS(export,file=file.path("results",output_name))
      }
    }
  }
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  run_highd(args[1])
}