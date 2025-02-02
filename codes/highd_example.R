rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(dplyr)
library(readr)

##### import functions #####
Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp")
source("functions/r_functions.R")


##### Read file for parameters #####
parameters <- as.data.frame(read_csv("results/simulation_details_highd.csv"))

#### Prompt to choose which simulation to run
# writeLines("You can write various IDs separated by commas")
id_chosen <- as.numeric(readline('Choose id:'))

sim_chosen <- parameters |> filter(id==id_chosen)
if(nrow(sim_chosen)!=1){print(paste0("Error: id ",id_chosen," doesn't exist or there's more than one")); next;}
# Parameters for all algorithms
total_simulations <- sim_chosen$simulations
temperatures <- as.numeric(sim_chosen[paste0("t",1:3)])
bal_f <- as.character(sim_chosen[paste0("bf",1:3)])
defined_seed <- sim_chosen$seed
set.seed(defined_seed)
#Parameters for PT with IIT
total_iter <- sim_chosen$iterations #300000 #Total number of steps to perform in each replica
iterswap <- sim_chosen$interswap #Total iterations before trying a replica swap
#Parameters for PT with a-IIT
sample_inter_swap <- sim_chosen$interswap #Number of original samples to get before trying a replica swap
total_swap <- sim_chosen$total_swap #Total number of swaps to try

# start_state <- sim_chosen$start_state;
alg <- sim_chosen$algorithm

export <- list();
#### Function depending on algorithm to use

writeLines(c("Parameters:",paste0("Algorithm: ",alg),
             paste0("Seed: ",defined_seed),
             paste0("Total simulations: ",total_simulations),
             paste0("Temperatures: ",paste(temperatures,collapse=',')),
             paste0("Balancing functions: ",paste(bal_f,collapse = ',')),
             paste0("Total iterations: ",total_iter),
             paste0("Try swaps:",iterswap),
             paste0("Samples in-between swaps: ",sample_inter_swap),
             paste0("Total swaps:",total_swap),"","","Confirm below"))

# check <- as.numeric(readline('ok? 1 Yes/ 0 No'))
check <- 1;
if(check!=1){print("modify parameters")}else{
  if(alg=="IIT"){
    # Only IIT
    # output_name <- paste0("IIT_","sim_",total_simulations,"_iter_",total_iter,"_s_",defined_seed,".Rds");
    output <- PT_IIT_sim(p,startsim=1, endsim=total_simulations,numiter=total_iter,iterswap=total_iter+1,temp=temperatures[1],bal_function=bal_f[1], bias_fix = TRUE,initial_state = start_state)
  }else{
    
    if(alg=="PT_IIT_Z"){
      # output_name <- paste0("PT_IIT_Z_","sim_",total_simulations,"_iter_",total_iter,"_iterswap_",iterswap,"_s_",defined_seed,".Rds");
      # Using Z factor bias correction
      output <- PT_IIT_sim(p,startsim=1, endsim=total_simulations,numiter=total_iter,iterswap,temperatures,bal_f, bias_fix = TRUE,initial_state = start_state)
      #round trip rate (NA for IIT)
      export[["round_trips"]] <- PT_RT(output[["ip"]], floor(total_iter/iterswap),total_simulations)
    }
    if(alg=="PT_IIT_no_Z"){
      # output_name <- paste0("PT_IIT_no_Z_","sim_",total_simulations,"_iter_",total_iter,"_iterswap_",iterswap,"_s_",defined_seed,".Rds");
      # Without Z factor bias correction
      output <- PT_IIT_sim(p,startsim=1, endsim=total_simulations,numiter=total_iter,iterswap,temperatures,bal_f, bias_fix = FALSE,initial_state = start_state)
      #round trip rate (NA for IIT)
      export[["round_trips"]] <- PT_RT(output[["ip"]], floor(total_iter/iterswap),total_simulations)
    }
    if(alg=="PT_A_IIT"){
      # Using A-IIT in each replica
      # output_name <- paste0("PT_A_IIT_","sim_",total_simulations,"_interswap_",sample_inter_swap,"_totalswap_",total_swap,"_s_",defined_seed,".Rds");
      output <- PT_a_IIT_sim(p,startsim=1, endsim=total_simulations,total_swap,sample_inter_swap,temperatures,bal_f,initial_state = start_state)
      #Number of iterations needed between swaps for each replica
      export[["total_iter"]] <- output[["total_iter"]]
      #round trip rate (NA for IIT)
      export[["round_trips"]] <- PT_RT(output[["ip"]],total_swap,total_simulations)
    }
    # Replica swap acceptance rate (NA for IIT)
    export[["swap_rate"]] <- output[["swap_rate"]]
    
  }
  
  export[["states"]] <- output[["states"]]
  export[["loglik_visited"]] <- output[["loglik_visited"]]
  export[["iter_visit"]]<- output[["iter_visit"]]
  output_name <- paste0("sim_highdim_id_",id_chosen,".Rds")
  saveRDS(export,file=file.path("results",output_name))
}