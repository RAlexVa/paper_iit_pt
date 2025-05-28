rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(dplyr, warn.conflicts = F)
library(readr)

find_temps <- function(list_ids){
  if(!("./results" %in% list.dirs(recursive=F))){## If there's no folder for exports
    stop("Wrong directory. There's no results folder for the output")
  }
  
  
    if(missing(list_ids)){## If parameters are not provided
      #### Prompt to choose which simulation to run
      writeLines("You can write various IDs separated by commas")
      list_ids <- readline('Choose id:')
    }
    if(!is.character(list_ids)){ list_ids <- as.character(list_ids)}
    list_ids <- as.numeric(unlist(strsplit(list_ids,",")))
    
    Rcpp::sourceCpp("functions/find_temp_func.cpp")# Source CPP functions
##### Read file for parameters #####
    parameters <- as.data.frame(read_csv("results/find_temps.csv", col_types = cols()))    

##### Start process for algorithms        
    for(id_chosen in list_ids){
      sim_chosen <- parameters |> filter(id==id_chosen)
      if(nrow(sim_chosen)!=1){print(paste0("Error: id ",id_chosen," doesn't exist or there's more than one")); next;}

      # Parameters for all algorithms
      alg <- sim_chosen$algorithm # algorithm
      interswap <- sim_chosen$interswap #Total iterations before trying a replica swap
      defined_seed <- sim_chosen$seed
      set.seed(defined_seed)
      p <- as.numeric(sim_chosen$p)
      theta <- as.numeric(sim_chosen$theta)
      bal_f <- as.character(sim_chosen$bf)
      model <- as.character(sim_chosen$model)
      temp_ini <- as.numeric(sim_chosen$t1)
      number_temperatures <- as.numeric(sim_chosen$num_temp)
      export <- list();
      #### Function depending on algorithm to use
#temperature_PT_IIT(int p,int interswap, double temp_ini, const std::string bal_function, const double& theta)      
      writeLines(c("Parameters:",
                   paste0("ID: ",id_chosen),
                   paste0("Algorithm: ",alg),
                   paste0("Problem: ",model),
                   paste0("Dimension: ",p),
                   paste0("Theta: ",theta),
                   paste0("Samples in-between swaps: ",interswap),
                   paste0("Initial temperature: ",temp_ini),
                   paste0("#Temperatures to find: ",number_temperatures),
                   paste0("Balancing function: ",bal_f)))
      
      
  for(t in 1:number_temperatures){
    if(alg=="PT_IIT_Z"){
      # Using Z factor bias correction
      #temperature_PT_IIT(int p,int interswap, double temp_ini, const std::string bal_function, const double& theta)
      output <- temperature_PT_IIT(p,interswap,temp_ini,bal_f, theta)
    }
    if(alg=="PT_A_IIT"){
      # Using A-IIT with multiplicity list in each replica
      # output_name <- paste0("PT_A_IIT_","sim_",total_simulations,"_interswap_",sample_inter_swap,"_totalswap_",total_swap,"_s_",defined_seed,".Rds");
      # output <- temperature_PT_a_IIT(p,1, total_simulations,total_swap,sample_inter_swap,burnin_iter,temperatures,bal_f,start_state,reduc_constant,reduc_model)
    }
    
    if(t==1){
      write(paste0("alg: ",alg,"\ntemp_ini: ",temp_ini,"\ntemp ",t+1," ",output), file = paste0("results/temperatures_id_ ",id_chosen,"_alg_",alg,".txt"), append = FALSE)
    }else{
      write(paste0("temp ",t+1," ",output), file = paste0("results/temperatures_id_ ",id_chosen,"_alg_",alg,".txt"), append = TRUE)
    }
    temp_ini <- output;
  } ## End for loop for the temperatures        
   
       
      
             
        
  }## End of the for loop that runs over IDS


  
}## End of function

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  run_lowd(args[1])
}




