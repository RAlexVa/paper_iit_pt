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
    
    
##### Read file for parameters #####
    parameters <- as.data.frame(read_csv("results/find_temps.csv", col_types = cols()))    
list_of_algs <- unique(parameters |> filter(id %in% list_ids) |> pull(algorithm))
# if(length(list_of_algs)>1){
#   stop("You cannot run more than 1 algorithm at a time");
# }else{
#   # if(list_of_algs=="PT_IIT_Z"){Rcpp::sourceCpp("functions/find_temp_parallel.cpp");}
#   # if(list_of_algs=="PT_A_IIT"){Rcpp::sourceCpp("functions/find_temp_func.cpp")}
# }
# Source CPP functions
Rcpp::sourceCpp("functions/find_temp_parallel.cpp");
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
      gibbs_step <- as.numeric(sim_chosen$gibbs_step)
      search_direction <- as.numeric(sim_chosen$direction)
      #### Function depending on algorithm to use
#temperature_PT_IIT(int p,int interswap, double temp_ini, const std::string bal_function, const double& theta)      
      writeLines(c("Parameters:",
                   paste0("ID: ",id_chosen),
                   paste0("Algorithm: ",alg),
                   paste0("Problem: ",model),
                   paste0("Dimension: ",p),
                   paste0("Theta: ",theta),
                   # paste0("Samples in-between swaps: ",interswap),
                   paste0("Gibbs Steps: ",gibbs_step),
                   paste0("Initial temperature: ",temp_ini),
                   paste0("#Temperatures to find: ",number_temperatures),
                   paste0("Direction: ",search_direction),
                   paste0("Balancing function: ",bal_f)))
      
      ##Transform balancing function to integer
      if(bal_f=="min"){bal_f=1;}else{bal_f=2;}
      
  for(t_counter in 1:number_temperatures){
    if(alg=="PT_IIT_Z"){
      # Using Z factor bias correction
      #find_temp_gibbs_PT_IIT(int p, int burn_in,double temp_ini, int bal_func, const double& theta, int gibbs_steps)
      output_list <- find_temp_gibbs_PT_IIT(p,burn_in=10,temp_ini,bal_f, theta,gibbs_step,search_direction)
      output <- output_list[["temp"]];
    }
    if(alg=="PT_A_IIT"){
      # Using A-IIT with multiplicity list in each replica
      #find_temp_gibbs_A_IIT(int p,int interswap, int burn_in,double temp_ini, int bal_func, const double& theta, int base_seed)
      output_list <- find_temp_gibbs_A_IIT(p,1,burn_in=10,temp_ini,bal_f,theta,defined_seed,search_direction)
      output <- output_list[["temp"]];
    }
    
    if(t_counter==1){
      write(paste0("alg: ",alg,"\ntemp_ini: ",temp_ini,
                   "\ntemp ",t_counter+1,": ",output,
                   ", swap_rate: ",output_list[["swap_rate"]],
                   ", num_swap: ",output_list[["swap"]],
                   ", seconds: ",output_list[["seconds"]]), 
            file = paste0("results/temperatures_id_",id_chosen,".txt"), 
            append = FALSE)
    }else{
      write(paste0("temp ",t_counter+1,": ",output,
                   ", swap_rate: ",output_list[["swap_rate"]],
                   ", num_swap: ",output_list[["swap"]],
                   ", seconds: ",output_list[["seconds"]]), 
            file = paste0("results/temperatures_id_",id_chosen,".txt"), 
            append = TRUE)
    }
## Output de number of iterations
    # if(alg=="PT_A_IIT"){
    #   write(paste0("Iterations: ",paste0(round(t(output_list[["iter"]]),3),collapse = ", ")),
    #         file = paste0("results/temperatures_id_",id_chosen,".txt"), 
    #         append = TRUE)}

    
    temp_ini <- output;
  } ## End for loop for the temperatures        
   
  }## End of the for loop that runs over IDS

}## End of function

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  find_temps(args[1])
}