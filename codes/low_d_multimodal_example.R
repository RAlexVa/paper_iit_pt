rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(dplyr)

##### import functions #####
Rcpp::sourceCpp("functions/cpp_functions.cpp")
source("functions/r_functions.R")


##### Define parameters #####

# Parameters for all algorithms
set.seed(153)
total_simulations <- 100#100
temperatures <- c(1,0.18,0.09,.001)
bal_f <- c("sq","sq","sq","sq")

#Parameters for PT with IIT
total_iter <- 500000 #300000 #Total number of steps to perform in each replica
iterswap <- 2000 #Total iterations before trying a replica swap

#Parameters for PT with a-IIT
sample_inter_swap <- 2000 #Number of original samples to get before trying a replica swap
total_swap <- 200 #Total number of swaps to try

##### Low-dimensional multimodal example #####

### Setup
p <- 16 #dimension
theta <- 8 #tail weight parameter

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
  pi.true[i+1] <-  ll_comp(NumberToVector(i,p),modes_list,theta,p)
}
pi.true <- exp(pi.true)/(length(modes_list)*(1+exp(-theta))^p)
# pi.true[modes+1] #True probability of the modes
# sum(pi.true[modes+1]) #Accumulated probability of the modes
# plot(pi.true)
# plot(pi.true[-(modes+1)])

#### Prompt to choose algorithm
writeLines(c('Choose algorithm:','1 is PT-IIT with Z factor correction','
2 is PT-IIT without Z factors (no bias correction)','
3 is PT with A-IIT in each replica','
4 is just IIT'))

alg <- as.numeric(readline('Select algortihm'))

export <- list();
#### Function depending on algorithm to use
if(alg %in% c(1,2,3,4)){
  writeLines(c("Parameters:",paste0("Total simulations: ",total_simulations),
               paste0("Temperatures: ",paste(temperatures,collapse=',')),
               paste0("Balancing functions: ",paste(bal_f,collapse = ',')),
               paste0("Total iterations: ",total_iter),
               paste0("Try swaps:",iterswap),
               paste0("Samples in-between swaps: ",sample_inter_swap),
               paste0("Total swaps:",total_swap),"","","Confirm below"))
  
  check <- as.numeric(readline('ok? 1 Yes/ 0 No'))
  if(check!=1){alg <- 0;print("modify parameters")}
if(alg==4){
  # Only IIT
  output_name <- paste0("IIT_","sim_",total_simulations,"_iter_",total_iter);
  output <- PT_IIT_sim(p,startsim=1, endsim=total_simulations,numiter=total_iter,iterswap=total_iter+1,temp=temperatures[1],bal_function=bal_f[1], bias_fix = TRUE)
}else{
  
if(alg==1){
  output_name <- paste0("PT_IIT_Z_","sim_",total_simulations,"_iter_",total_iter,"_iterswap_",iterswap);
  # Using Z factor bias correction
  output <- PT_IIT_sim(p,startsim=1, endsim=total_simulations,numiter=total_iter,iterswap,temperatures,bal_f, bias_fix = TRUE)
  #round trip rate (NA for IIT)
  export[["round_trips"]] <- PT_RT(output[["ip"]], floor(total_iter/iterswap),total_simulations)
}
if(alg==2){
  output_name <- paste0("PT_IIT_no_Z_","sim_",total_simulations,"_iter_",total_iter,"_iterswap_",iterswap);
  # Without Z factor bias correction
  output <- PT_IIT_sim(p,startsim=1, endsim=total_simulations,numiter=total_iter,iterswap,temperatures,bal_f, bias_fix = FALSE)
  #round trip rate (NA for IIT)
  export[["round_trips"]] <- PT_RT(output[["ip"]], floor(total_iter/iterswap),total_simulations)
}
if(alg==3){
  # Using A-IIT in each replica
  output_name <- paste0("PT_A_IIT_","sim_",total_simulations,"_interswap_",sample_inter_swap,"_totalswap_",total_swap);
  output <- PT_a_IIT_sim(p,startsim=1, endsim=total_simulations,total_swap,sample_inter_swap,temperatures,bal_f)
  #Number of iterations needed between swaps for each replica
  export[["total_iter"]] <- output[["total_iter"]]
  #round trip rate (NA for IIT)
  export[["round_trips"]] <- PT_RT(output[["ip"]],total_swap,total_simulations)
}
# Replica swap acceptance rate (NA for IIT)
export[["swap_rate"]] <- output[["swap_rate"]]

}

#Compute estimated density
output[["est_pi"]] <- t(t(output[["est_pi"]])/colSums(output[["est_pi"]])) 
# Total Variation Distance
export[["tvd"]] <- apply(output[["est_pi"]], 2,TVD,pi.est=pi.true)

# Time of first visit
export[["mode_visit"]] <- t(output[["visits"]][modes+1,])


save(export,file=file.path("results",output_name))
}else{ print("Choose a valid option")}







