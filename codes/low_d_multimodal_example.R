rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(dplyr)

##### import functions #####
Rcpp::sourceCpp("functions/cpp_functions.cpp")
source("functions/r_functions.R")

##### Low-dimensional multimodal example #####

##### simulations for multimodal example #####
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


pi.true <- rep(0,2^p)
for(i in 0:(2^p -1)){
  pi.true[i+1] <-  ll_comp(NumberToVector(i,p),modes_list,theta,p)
}
#True probability
pi.true <- exp(pi.true)/(length(modes_list)*(1+exp(-theta))^p)

pi.true[modes+1] #True probability of the modes
sum(pi.true[modes+1]) #Accumulated probability of the modes
# plot(pi.true)
# plot(pi.true[-(modes+1)])

##### IIT with PT simulations #####
Rcpp::sourceCpp("functions/cpp_functions.cpp")
total_iter <- 300000
iterswap <- 2000
total_simulations <- 10
temperatures <- c(1,0.18,0.09,.001)
bal_f <- c("sq","sq","sq","sq")
set.seed(153)
ex1 <- PT_IIT_sim(p,startsim=1, endsim=total_simulations,numiter=total_iter,iterswap,temperatures,bal_f, bias_fix = TRUE)
ex2 <- PT_IIT_sim(p,startsim=1, endsim=total_simulations,numiter=total_iter,iterswap,temperatures,bal_f, bias_fix = FALSE)
# Only IIT
ex3 <- PT_IIT_sim(p,startsim=1, endsim=total_simulations,numiter=total_iter,iterswap=total_iter+1,temp=c(1),bal_function=c("sq"), bias_fix = TRUE)

output <- ex1

#Compute estimated density
output[["est_pi"]] <- t(t(output[["est_pi"]])/colSums(output[["est_pi"]])) 
# Total Variation Distance
apply(output[["est_pi"]], 2,TVD,pi.est=pi.true)

# Time of first visit
output[["visits"]][modes+1,]

# Replica swap acceptance rate
output[["swap_rate"]]

#round trip rate
PT_RT(output[["ip"]], total_iter,iterswap,total_simulations)



##### Adaptive PT simulations #####
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("functions/cpp_functions.cpp")
sample_inter_swap <- 3000
total_swap <- 200
temperatures <- c(1,0.18,0.09,.001)
bal_f <- c("sq","sq","sq","sq")
set.seed(153)
ex1 <- PT_a_IIT_sim(p,startsim=1, endsim=2,total_swap,sample_inter_swap,temperatures,bal_f)

pi.est <- ex1[[1]]/sum(ex1[[1]])
TVD(pi.est,pi.true)
ex1[[3]][modes+1]
ex1[[4]]
ex1[[5]]
