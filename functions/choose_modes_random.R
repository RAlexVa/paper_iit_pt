rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(dplyr, warn.conflicts = F)
library(readr)

##### import functions #####

source("functions/r_functions.R");
Rcpp::sourceCpp("functions/highdim_2_parallel.cpp");



choose_modes <- function(p,num_temps,num_modes){
  Q_matrix <- create_mode_matrix(p,num_modes)
  
  dummy_vec <- Random_modes(p,num_temps,Q_matrix)
  
  # return(length(unique(dummy_vec)))
  return(dummy_vec)
}


### Start running the code


#For 2 modes we use 7 replicas
#For 5 modes we use 13 replicas
#For 7 modes we use 25 replicas

start <- 1
end <- 500
p <- 1000
vec_temps <- c(7,13,25)
vec_modes <- c(2,5,7)

storage <- matrix(0,nrow=end-start+1,ncol=3)


for(m in 1:3){
  num_temps <- vec_temps[m]
  num_modes <- vec_modes[m]
  for(i in start:end){
    set.seed(i)
    temp <- choose_modes(p,num_temps,num_modes)
    storage[i,m] <- length(unique(temp))
  }
}


checking_5_modes <- c(104,105,110,114,116,119,120,121,122,126,131,137,138,140,141,146,152,154,155,156,157,158,159,160,161,168,170,173,175,177,180,181,184,185,186,187,189,190,192,193,194,196,200,203,206,207,208,210,212,215,217,219,220,226,227,228,229,230,232,235,237,240,243,244,247,249)

storage[sim_not_visiting,]
#No está replicando las modas que escoje en las simulaciones, 
#quizá los cáclulos que hace antes interfieren o no se
#El punto es que no pitufó

