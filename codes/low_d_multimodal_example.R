rm(list=ls())
library(Rcpp)
library(RcppArmadillo)

##### Some functions #####
NumberToVector <- function(number,N)
{
  if(number < 0 || number >= 2^N)
  {
    print("Error")
    return(0)
  }
  X = rep(0, N)
  # number = number - 1
  for(k in 1 : N)
  {
    X[k] = number %% 2
    number = as.integer(number / 2)
  }
  return(X)
}
#Replica swap
PT_swap_rate <- function(index_process){
  
  head(index_process,10)
  ch1 <- index_process[c(T,F),]
  swap_rate1 <- colSums(index_process[-1, ] != index_process[-nrow(ip1), ])
  
  #Need to multiply by 2 since not every swap is trying to swap all temperatures
  return(2*swap_rate1/(nrow(index_process)-1));
}
#Total variation distance
TVD <- function(pi,pi.est){
  return(0.5*sum(abs(pi-pi.est)))
}
##### Low-dimensional multimodal example #####
Rcpp::sourceCpp("functions/cpp_functions.cpp")

### computing normalizing constant

### modes definition
d=8
mod1 <- rep(1,d)
mod2 <- c(1,0,1,0,1,0,1,0)
mod3 <- c(0,1,0,1,0,1,0,1)
mod4 <- c(1,1,1,1,0,0,0,0)
mod5 <- c(0,0,0,0,1,1,1,1)
mod6 <- c(1,0,0,0,0,0,0,1)
mod7 <- c(0,0,0,1,1,0,0,0)

# ll_comp_2 <- function(X,theta=1){
#   if(length(X)!=d){print('Error in dimension');return(-1);}
#   total=0;
#   for(m in 1:7){# loop for modes
#     total = total+theta*sum(abs(X-eval(parse(text=paste0("mod",m)))))
#   }
# return(-total);
# }


# Testing function
NumberToVector(2^d-1,d)
NumberToVector(0,d)
NumberToVector(7,d)
##### corresponding states for the modes
# vec_to_num(mod1)
# vec_to_num(mod2)
# vec_to_num(mod3)
# vec_to_num(mod4)
# vec_to_num(mod5)
# vec_to_num(mod6)
# vec_to_num(mod7)
modes <- c(vec_to_num(mod1),  vec_to_num(mod2),  vec_to_num(mod3),  vec_to_num(mod4),  vec_to_num(mod5),  vec_to_num(mod6),  vec_to_num(mod7))
#Computing log-likelihoods for all states
loglik_vec <- rep(0,2^d)
for(i in 0:(2^d -1)){
  loglik_vec[i+1] <-  ll_comp(NumberToVector(i,d))
}

#Confirming the log-likelihood of the probabilities
loglik_vec[modes+1]
min(loglik_vec[modes+1])

#Confirming those are the only modes
modes[modes>=min(loglik_vec[modes+1])]

#computing normalizing constant and comparing with theoretical constant
sum(exp(loglik_vec))
7*(1+exp(-1))^d


##### Testing loglik function from CPP #####
loglik(mod1)
ll_comp(mod1)
loglik(c(1,1,1,0,0,0,0,0))
ll_comp(c(1,1,1,0,0,0,0,0))

##### simulations for multimodal example #####
#For p=2
p <- 16 #dimension
theta <- 15
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
ll_comp <- function(X, theta=1){
  if(length(X)!=p){print('Error in dimension');return(-1);}
  total=0;
  for(m in 1:length(modes_list)){# loop for modes
    total = total+exp(-theta*sum(abs(X-modes_list[[m]])))
  }
  return(log(total));
}
pi.true <- rep(0,2^p)
for(i in 0:(2^p -1)){
  pi.true[i+1] <-  ll_comp(NumberToVector(i,p),theta)
}
#True probability
pi.true <- exp(pi.true)/(length(modes_list)*(1+exp(-theta))^p)

pi.true[modes+1]
sum(pi.true[modes+1])


##### IIT with PT simulations #####
Rcpp::sourceCpp("functions/cpp_functions.cpp")
total_iter <- 300000
iterswap <- 2000
temperatures <- c(1,0.18,0.09,.001)
bal_f <- c("sq","sq","sq","sq")
set.seed(153)
ex1 <- PT_IIT_sim(p,startsim=1, endsim=1,numiter=total_iter,iterswap,temperatures,bal_f, bias_fix = TRUE)
ex2 <- PT_IIT_sim(p,startsim=1, endsim=1,numiter=total_iter,iterswap,temperatures,bal_f, bias_fix = FALSE)
ex3 <- PT_IIT_sim(p,startsim=1, endsim=1,numiter=total_iter,iterswap=1000000,temp=c(1),bal_function=c("sq"), bias_fix = TRUE)

pi.est <- ex1[[1]]/sum(ex1[[1]])
TVD(pi.est,pi.true)
ex1[[4]]
ex1[[3]][modes+1]
pi.est[modes+1]

pi.est <- ex2[[1]]/sum(ex2[[1]])
TVD(pi.est,pi.true)
ex2[[4]]
ex2[[3]][modes+1]
pi.est[modes+1]


pi.est <- ex3[[1]]/sum(ex3[[1]])
TVD(pi.est,pi.true)
ex3[[3]][modes+1]
pi.est[modes+1]

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
