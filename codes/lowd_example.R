rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(dplyr)
library(readr)

##### import functions #####
# Rcpp::sourceCpp("functions/cpp_functions.cpp")
source("functions/r_functions.R")
run_lowd <- function(list_ids){
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
    parameters <- as.data.frame(read_csv("results/simulation_details_lowd.csv"))
    {#Code to compute the true pi for both models and have it in memory
      Rcpp::sourceCpp("functions/cpp_functions.cpp")# To use vec_to_num function
      ##### Low-dimensional 7 modes setup #####
      {
        p <- 16 #dimension
        theta <- 15 #tail weight parameter
        
        # Modes definition
        mod1 <- rep(1,p)
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
        pi.true.7 <- pi.true;
        modes.7 <- modes;
      }
      ##### Low-dimensional BI-modal setup #####
      {
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
        pi.true.2 <- pi.true;
        modes.2 <- modes;
      }
      
      
    }
    #Check how many models we're doing
    tot_models <- unique(parameters|> filter(id %in% list_ids) |> pull(model))
    if(length(tot_models)==1){only_1_model <- TRUE;}else{only_1_model <- FALSE}
    # In case only 1 model is chosen, we only read 1 model for all the IDs    
    if(only_1_model){
      print("Reading one set of C++ functions")
      if(tot_models=="7_mode"){
        Rcpp::sourceCpp("functions/cpp_functions.cpp")
        pi.true <- pi.true.7;
        modes <- modes.7;
        }
      if(tot_models=="bimodal"){
        Rcpp::sourceCpp("functions/cpp_functions_lowdim_2.cpp")
        pi.true <- pi.true.2;
        modes <- modes.2;
        }
    }
    #Start process for algorithms        
    for(id_chosen in list_ids){
      sim_chosen <- parameters |> filter(id==id_chosen)
      if(nrow(sim_chosen)!=1){print(paste0("Error: id ",id_chosen," doesn't exist or there's more than one")); next;}
      #In case there are more than 1 model, I need to re-read functions depending on model
      if(!only_1_model){
        print(paste0("Reading C++ functions for id: ",id_chosen," model: ",sim_chosen$model))
        if(tot_models=="7_mode"){
          Rcpp::sourceCpp("functions/cpp_functions.cpp")
          pi.true <- pi.true.7;
          modes <- modes.7;
        }
        if(tot_models=="bimodal"){
          Rcpp::sourceCpp("functions/cpp_func_multilow.cpp")
          pi.true <- pi.true.2;
          modes <- modes.2;
        }
      }      
      # Parameters for all algorithms
      total_simulations <- sim_chosen$simulations
      temperatures <- as.numeric(sim_chosen |> select(matches("^t\\d+$")))
      temperatures <- temperatures[!is.na(temperatures)] #Ignore NA temperatures
      bal_f <- as.character(sim_chosen|> select(matches("^bf")))
      bal_f <- bal_f[!is.na(bal_f)] #Ignore NA balancing functions
      if(length(bal_f)==0){bal_f <- "sq"}
      burnin_iter <- sim_chosen$burn_in #Number of iterations for burn-in
      start_state <- sim_chosen$start_state;#starting point
      alg <- sim_chosen$algorithm # algorithm
      defined_seed <- sim_chosen$seed
      set.seed(defined_seed)
      #Parameters for PT with IIT
      total_iter <- sim_chosen$iterations #300000 #Total number of steps to perform in each replica
      iterswap <- sim_chosen$interswap #Total iterations before trying a replica swap
      #Parameters for PT with a-IIT
      sample_inter_swap <- sim_chosen$interswap #Number of original samples to get before trying a replica swap
      total_swap <- sim_chosen$total_swap #Total number of swaps to try
      reduc_constant <- sim_chosen$reduc_constant
      reduc_model <- sim_chosen$reduc_model
      
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
          #PT_IIT_sim(int p,int startsim,int endsim, int numiter,int iterswap,int burn_in, vec temp, const std::vector<std::string>& bal_function, bool bias_fix, int initial_state)
          output <- PT_IIT_sim(p,startsim=1, endsim=total_simulations,numiter=total_iter,iterswap=total_iter+1,burnin_iter,temp=temperatures[1],bal_function=bal_f[1], bias_fix = TRUE,initial_state = start_state)
        }else{
          
          if(alg=="PT_IIT_Z"){
            # Using Z factor bias correction
            #PT_IIT_sim(int p,int startsim,int endsim, int numiter,int iterswap,int burn_in, vec temp, const std::vector<std::string>& bal_function, bool bias_fix, int initial_state)
            output <- PT_IIT_sim(p,1, total_simulations,total_iter,iterswap,burnin_iter,temperatures,bal_f, TRUE,start_state)
            #round trip rate (NA for IIT)
            export[["round_trips"]] <- PT_RT(output[["ip"]], floor(total_iter/iterswap),total_simulations)
          }
          if(alg=="PT_IIT_no_Z"){
            # output_name <- paste0("PT_IIT_no_Z_","sim_",total_simulations,"_iter_",total_iter,"_iterswap_",iterswap,"_s_",defined_seed,".Rds");
            # Without Z factor bias correction
            output <- PT_IIT_sim(p,1, total_simulations,total_iter,iterswap,burnin_iter,temperatures,bal_f, FALSE,start_state)
            #round trip rate (NA for IIT)
            export[["round_trips"]] <- PT_RT(output[["ip"]], floor(total_iter/iterswap),total_simulations)
          }
          if(alg=="PT_A_IIT"){
            # Using A-IIT with multiplicity list in each replica
            # output_name <- paste0("PT_A_IIT_","sim_",total_simulations,"_interswap_",sample_inter_swap,"_totalswap_",total_swap,"_s_",defined_seed,".Rds");
            output <- PT_a_IIT_sim(p,1, total_simulations,total_swap,sample_inter_swap,burnin_iter,temperatures,bal_f,start_state,reduc_constant,reduc_model)
            #Number of iterations needed between swaps for each replica
            export[["total_iter"]] <- output[["total_iter"]]
            #round trip rate (NA for IIT)
            export[["round_trips"]] <- PT_RT(output[["ip"]],total_swap,total_simulations)
          }
          if(alg=="PT_A_IIT_RF"){
            # Using A-IIT with weights in each replica
            #PT_a_IIT_sim_RF(int p,int startsim,int endsim, int numiter,int iterswap,int burn_in, vec temp, const std::vector<std::string>& bal_function, bool bias_fix, int initial_state)
            output <- PT_a_IIT_sim_RF(p,1,total_simulations,total_iter,iterswap,burnin_iter,temperatures,bal_f,TRUE,start_state,reduc_constant,reduc_model)
            #round trip rate (NA for IIT)
            export[["round_trips"]] <- PT_RT(output[["ip"]], floor(total_iter/iterswap),total_simulations)
          }
          # Replica swap acceptance rate (NA for IIT)
          export[["swap_rate"]] <- output[["swap_rate"]]
        }
        
        #Compute estimated density
        export[["est_pi"]] <- t(t(output[["est_pi"]])/colSums(output[["est_pi"]]))
        # export[["est_pi"]] <- output[["est_pi"]]
        # Total Variation Distance
        export[["tvd"]] <- apply(output[["est_pi"]], 2,TVD,pi.est=pi.true)
        # estimated density on the modes
        export[["pi_modes"]] <- output[["est_pi"]][modes+1,]
        # Time of first visit
        export[["mode_visit"]] <- t(output[["visits"]][modes+1,])
        #Running time
        export[["time_taken"]] <- output[["time_taken"]]
        #Index process
        export[["ip"]] <- output[["ip"]]
        output_name <- paste0("sim_lowdim_id_",id_chosen,".Rds")
        saveRDS(export,file=file.path("results",output_name))
      }
    }
  }
}





