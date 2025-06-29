rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(dplyr, warn.conflicts = F)
library(readr)

##### import functions #####

source("functions/r_functions.R");
Rcpp::sourceCpp("functions/highdim_2_parallel.cpp");
find_temp_highd <- function(list_ids){
  target_sr <- 0.234
  
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
    parameters <- as.data.frame(read_csv("inputs/find_temps.csv", col_types = cols()))
    
    

    #Start process for algorithms
    for(id_chosen in list_ids){
      sim_chosen <- parameters |> filter(id==id_chosen)
      if(nrow(sim_chosen)!=1){print(paste0("Error: id ",id_chosen," doesn't exist or there's more than one")); next;}
      #In case there are more than 1 model, I need to re-read functions depending on model
 
      # Parameters for all algorithms
      total_simulations <- sim_chosen$tot_sim
      p <- sim_chosen$p
      num_temps <- sim_chosen$num_temp
      num_modes <- sim_chosen$num_modes
      search_direction <- sim_chosen$direction
      temp_ini<- sim_chosen$t1
      temperatures <- rep(0,num_temps)
      temperatures[1] <- temp_ini
      initial_rho <- sim_chosen$rho
      rho_vector <- rep(initial_rho,num_temps-1)
      for(t in 2:num_temps){
        if(sign(search_direction)<0){temperatures[t]=temperatures[t-1]/(1+exp(rho_vector[t-1]))}
        if(sign(search_direction)>0){temperatures[t]=temperatures[t-1]*(1+exp(rho_vector[t-1]))}
        
      }
      original_temperatures <- temperatures;
      bal_f <- as.character(sim_chosen$bf)
      burnin_iter <- sim_chosen$burn_in #Number of iterations for burn-in
      alg <- sim_chosen$algorithm
      defined_seed <- sim_chosen$seed
      set.seed(defined_seed)
      #Parameters for PT with IIT
      total_iter <- sim_chosen$iterations 
      #Parameters for PT with A-IIT
      sample_inter_swap <- sim_chosen$interswap #Number of original samples to get before trying a replica swap
      total_swap <- sim_chosen$total_swap #Total number of swaps to try
     theta_chosen <- sim_chosen$theta;
      export <- list();
      #### Function depending on algorithm to use
      
      writeLines(c("Parameters:",
                   paste0("ID: ",id_chosen),
                   paste0("Algorithm: ",alg),
                   paste0("Problem: ",sim_chosen$model),
                   paste0("Num Modes: ",num_modes),
                   paste0("Number of Replicas: ",length(temperatures)),
                   paste0("Theta: ",theta_chosen),
                   paste0("Total simulations: ",total_simulations),
                   paste0("Burn-in iterations: ",burnin_iter),
                   paste0("Total iterations: ",total_iter),
                   paste0("Temperatures: ",paste(temperatures,collapse=',')),
                   paste0("Search direction: ",search_direction),
                   paste0("Balancing function: ",paste(bal_f,collapse = ',')),
                   paste0("Try swaps:",sample_inter_swap),
                   paste0("Samples in-between swaps: ",sample_inter_swap),
                   paste0("Total swaps:",total_swap)))
      
      #re-create the balancing function vector
      #All replicas use the same balancing function
      # bal_f <- rep(bal_f,length(temperatures))
      if(bal_f=="min"){bal_f=1;}else{bal_f=2;}

  check_convergence <- rep(F,length(rho_vector))

  check_loop <- TRUE;
  while(check_loop){
    if(alg=="PT_IIT_Z"){
      # Using Z factor bias correction
      # PT_IIT_sim(int p,int startsim,int endsim, int numiter, int iterswap,int burn_in, vec temp, int bal_func, bool bias_fix,const std::string& filename,int num_states_visited,const std::vector<int>& starting_coord, double theta)
      output <- PT_IIT_sim(p,1,total_simulations,total_iter,sample_inter_swap,burnin_iter,temperatures,bal_f,TRUE,"",0,0,theta_chosen,num_modes)
    }
    if(alg=="PT_A_IIT"){
      # Using A-IIT in each replica
      # PT_a_IIT_sim(int p,int startsim,int endsim, int total_swaps,int sample_inter_swap,int burn_in, vec temp, const int bal_func,const std::string& filename,int num_states_visited,const std::vector<int>& starting_coord, double decreasing_constant,std::string reduc_model, double theta)
      output <- PT_a_IIT_sim(p,1,total_simulations,total_swap,sample_inter_swap,burnin_iter,temperatures,bal_f,"",0,0,0,"never",theta_chosen,num_modes)
    }

    
    swap_rate <- as.data.frame(output[["swap_rate"]])
    writeLines(c("swap rate: "));
    print(swap_rate);
    summary_sr <- colSums(swap_rate)/nrow(swap_rate) #Get avg. swap rate
    distance_sr <- summary_sr-target_sr;
    # new_rho <- rho_vector+3*distance_sr;
    new_rho <- rho_vector+distance_sr;
    for(i in 1:length(summary_sr)){
      if(summary_sr[i]>0.234 && summary_sr[i]<0.27){
        #Dont change the rho_factor
        check_convergence[i] <- TRUE
      }else{
        rho_vector[i] <- new_rho[i];
      }
    }
    #Update temperatures
    for(t in 2:num_temps){
      if(sign(search_direction)<0){temperatures[t]=temperatures[t-1]/(1+exp(rho_vector[t-1]))}
      if(sign(search_direction)>0){temperatures[t]=temperatures[t-1]*(1+exp(rho_vector[t-1]))}
      
    }
    
 print(paste0("Rho Vector ",paste0(rho_vector,collapse=",")));
    

    writeLines(c(paste0("ID: ",id_chosen),
                 paste0("Temperatures: ",paste(round(temperatures,6),collapse=',')),
                 paste0("Swap Rate: ",paste(summary_sr,collapse = ',')),
                 paste0("Time: ",mean(output[["time_taken"]]))))
    
    if(all(check_convergence)){
      check_loop <- FALSE;
    }else{
      check_convergence <- rep(F,length(swap_rate))
    }
  }## End of while loop
  
  ## Export of results
  if(alg=="PT_IIT_Z"){
    write(paste0("alg: ",alg,"\ntemp_ini: ",temp_ini,
                 "\n temps: ",paste0(temperatures,collapse=","),
                 ",\n swap_rate: ",paste0(summary_sr,collapse=","),
                 ",\n seconds: ",mean(output[["time_taken"]])), 
          file = paste0("results/temperatures_id_",id_chosen,".txt"), 
          append = FALSE)
  }
  
  
  if(alg=="PT_A_IIT"){
    avg_iter <- rowSums(colSums(output[["total_iter"]]))/(total_swap*total_simulations)
    write(paste0("alg: ",alg,
                 "\ntemps: ",paste0(temperatures,collapse=","),
                 "\niterations: ",paste0(avg_iter,collapse=","),
                 ",\nswap_rate: ",paste0(summary_sr,collapse=","),
                 ",\nseconds: ",mean(output[["time_taken"]])), 
          file = paste0("results/temperatures_id_",id_chosen,".txt"), 
          append = FALSE)
  }
  
    }##End for loop for each ID chosen
  }## End else if there's no results folder
  
}### End of function

find_temp_highd_recurrent <- function(list_ids,min_temp_to_find=20){
  target_sr <- 0.234
  
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
    parameters <- as.data.frame(read_csv("inputs/find_temps.csv", col_types = cols()))
    
    
    havent_reach_min <- TRUE;
    after_first_run <- FALSE
    #Start process for algorithms
    for(id_chosen in list_ids){
      
      while(havent_reach_min){
      
      sim_chosen <- parameters |> filter(id==id_chosen)
      if(nrow(sim_chosen)!=1){print(paste0("Error: id ",id_chosen," doesn't exist or there's more than one")); next;}
      #In case there are more than 1 model, I need to re-read functions depending on model
      
      # Parameters for all algorithms
      total_simulations <- sim_chosen$tot_sim
      p <- sim_chosen$p
      num_temps <- sim_chosen$num_temp
      num_modes <- sim_chosen$num_modes
      search_direction <- sim_chosen$direction
      if(after_first_run){
        temp_ini <- new_temp_ini;
      }else{temp_ini<- sim_chosen$t1}

      temperatures <- rep(0,num_temps)
      temperatures[1] <- temp_ini
      initial_rho <- sim_chosen$rho
      rho_vector <- rep(initial_rho,num_temps-1)
      for(t in 2:num_temps){
        if(sign(search_direction)<0){temperatures[t]=temperatures[t-1]/(1+exp(rho_vector[t-1]))}
        if(sign(search_direction)>0){temperatures[t]=temperatures[t-1]*(1+exp(rho_vector[t-1]))}
        
      }
      original_temperatures <- temperatures;
      bal_f <- as.character(sim_chosen$bf)
      burnin_iter <- sim_chosen$burn_in #Number of iterations for burn-in
      alg <- sim_chosen$algorithm
      defined_seed <- sim_chosen$seed
      set.seed(defined_seed)
      #Parameters for PT with IIT
      total_iter <- sim_chosen$iterations 
      #Parameters for PT with A-IIT
      sample_inter_swap <- sim_chosen$interswap #Number of original samples to get before trying a replica swap
      total_swap <- sim_chosen$total_swap #Total number of swaps to try
      theta_chosen <- sim_chosen$theta;
      export <- list();
      #### Function depending on algorithm to use
      
      writeLines(c("Parameters:",
                   paste0("ID: ",id_chosen),
                   paste0("Algorithm: ",alg),
                   paste0("Problem: ",sim_chosen$model),
                   paste0("Num Modes: ",num_modes),
                   paste0("Number of Replicas: ",length(temperatures)),
                   paste0("Theta: ",theta_chosen),
                   paste0("Total simulations: ",total_simulations),
                   paste0("Burn-in iterations: ",burnin_iter),
                   paste0("Total iterations: ",total_iter),
                   paste0("Temperatures: ",paste(temperatures,collapse=',')),
                   paste0("Search direction: ",search_direction),
                   paste0("Balancing function: ",paste(bal_f,collapse = ',')),
                   paste0("Try swaps:",sample_inter_swap),
                   paste0("Samples in-between swaps: ",sample_inter_swap),
                   paste0("Total swaps:",total_swap)))
      
      #re-create the balancing function vector
      #All replicas use the same balancing function
      # bal_f <- rep(bal_f,length(temperatures))
      if(bal_f=="min"){bal_f=1;}else{bal_f=2;}
      
      check_convergence <- rep(F,length(rho_vector))
      
      check_loop <- TRUE;
      while(check_loop){
        if(alg=="PT_IIT_Z"){
          # Using Z factor bias correction
          # PT_IIT_sim(int p,int startsim,int endsim, int numiter, int iterswap,int burn_in, vec temp, int bal_func, bool bias_fix,const std::string& filename,int num_states_visited,const std::vector<int>& starting_coord, double theta)
          output <- PT_IIT_sim(p,1,total_simulations,total_iter,sample_inter_swap,burnin_iter,temperatures,bal_f,TRUE,"",0,0,theta_chosen,num_modes)
        }
        if(alg=="PT_A_IIT"){
          # Using A-IIT in each replica
          # PT_a_IIT_sim(int p,int startsim,int endsim, int total_swaps,int sample_inter_swap,int burn_in, vec temp, const int bal_func,const std::string& filename,int num_states_visited,const std::vector<int>& starting_coord, double decreasing_constant,std::string reduc_model, double theta)
          output <- PT_a_IIT_sim(p,1,total_simulations,total_swap,sample_inter_swap,burnin_iter,temperatures,bal_f,"",0,0,0,"never",theta_chosen,num_modes)
        }
        
        
        swap_rate <- as.data.frame(output[["swap_rate"]])
        writeLines(c("swap rate: "));
        print(swap_rate);
        summary_sr <- colSums(swap_rate)/nrow(swap_rate) #Get avg. swap rate
        distance_sr <- summary_sr-target_sr;
        # new_rho <- rho_vector+3*distance_sr;
        new_rho <- rho_vector+distance_sr;
        for(i in 1:length(summary_sr)){
          if(summary_sr[i]>0.234 && summary_sr[i]<0.27){
            #Dont change the rho_factor
            check_convergence[i] <- TRUE
          }else{
            rho_vector[i] <- new_rho[i];
          }
          # check_convergence[i] <- TRUE
        }
        #Update temperatures
        for(t in 2:num_temps){
          if(sign(search_direction)<0){temperatures[t]=temperatures[t-1]/(1+exp(rho_vector[t-1]))}
          if(sign(search_direction)>0){temperatures[t]=temperatures[t-1]*(1+exp(rho_vector[t-1]))}
          
        }
        
        print(paste0("Rho Vector ",paste0(rho_vector,collapse=",")));
        
        
        writeLines(c(paste0("ID: ",id_chosen),
                     paste0("Temperatures: ",paste(round(temperatures,6),collapse=',')),
                     paste0("Swap Rate: ",paste(summary_sr,collapse = ',')),
                     paste0("Time: ",mean(output[["time_taken"]]))))
        
        if(all(check_convergence)){
          check_loop <- FALSE;
        }else{
          check_convergence <- rep(F,length(swap_rate))
        }
      }## End of while loop
      
      ## Export of results
      if(alg=="PT_IIT_Z"){
        write(paste0("\nalg: ",alg,"\ntemp_ini: ",temp_ini,
                     "\n temps: ",paste0(temperatures,collapse=","),
                     ",\n swap_rate: ",paste0(summary_sr,collapse=","),
                     ",\n seconds: ",mean(output[["time_taken"]])), 
              file = paste0("results/temperatures_id_",id_chosen,".txt"), 
              append = TRUE)
      }
      
      
      if(alg=="PT_A_IIT"){
        avg_iter <- rowSums(colSums(output[["total_iter"]]))/(total_swap*total_simulations)
        write(paste0("\nalg: ",alg,
                     "\ntemps: ",paste0(temperatures,collapse=","),
                     "\niterations: ",paste0(avg_iter,collapse=","),
                     ",\nswap_rate: ",paste0(summary_sr,collapse=","),
                     ",\nseconds: ",mean(output[["time_taken"]])), 
              file = paste0("results/temperatures_id_",id_chosen,".txt"), 
              append = TRUE)
      }
      
      if(min(temperatures)<=min_temp_to_find){#We've reached the min_temperature
        havent_reach_min <- FALSE;
      }else{
        after_first_run <- TRUE
        new_temp_ini <- min(temperatures);
      }
    }# End while reach_min
    }##End for loop for each ID chosen
  }## End else if there's no results folder
  
}### End of function


if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  find_temp_highd(args[1])
}

