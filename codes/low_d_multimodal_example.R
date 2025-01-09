rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(dplyr)

##### Some functions #####
#Number to vector to compute TVD with target distribution
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
#Total variation distance
TVD <- function(pi,pi.est){
  if(sum(pi)!=1){pi <- pi/sum(pi)}
  if(sum(pi.est)!=1){pi.est <- pi.est/sum(pi.est)}
  return(0.5*sum(abs(pi-pi.est)))
}

#Compute roundtrip rate
PT_RT <- function(ip,total_iter,iterswap, total_sim){
  # ip <- ex2[["ip"]]+1
  ip <- as.data.frame(ip)
  total_swaps <- floor(total_iter/iterswap)
  if(nrow(ip)!=(total_sim*total_swaps+1)){print("Error with size of input")}
  
  max_temp <- max(ip[1,])
  min_temp <- min(ip[1,])
  
  initial_ip <- ip[1,]
  # current_state <- (initial_ip==max_temp) - (initial_ip==min_temp)
  # 
  # last_visited <- (initial_ip==max_temp | initial_ip==min_temp)*initial_ip
  # 
  # ip2 <- ip[-1,]
  ip <- ip |> slice(-1)
  ip$sim <- rep(1:total_sim,each=total_swaps)
  ip <- ip |> select(sim,everything())
  
#The round trip starts with the hottest temperature
  
  rt_count <- matrix(0,nrow = total_sim,ncol=ncol(ip)-1) # to track round trip rate
  
  for (c in 2:ncol(ip)){#Loop for each replica
    first_last <- (initial_ip[c-1]==max_temp)*max_temp + (initial_ip[c-1]==min_temp)*min_temp
    # last <- first_last
    for(i in 1:nrow(ip)){#Loop for each iteration
    # for(i in 1201:nrow(ip)){  
      if(i==1){# Restart the index process for each simulation
        last <- first_last
        start_rt <- FALSE
      }else if(ip[i,1]!=ip[i-1,1]){
        last <- first_last
        start_rt <- FALSE
      }
      current <- ip[i,c]#Current iteration
      if(current==max_temp & start_rt==FALSE){
        # print(paste0("Start round trip for simulation ",ip[i,1]," at iteration ", i));
        start_rt <- TRUE;
        last <- current;
            # if(current!=last){rt_count[ip[i,1],c-1]=rt_count[ip[i,1],c-1]-0.5;}
      next;
      }#We start counting RT from the max temp
      
      if(start_rt){
        if(current==min_temp|current==max_temp){ #If it's one of the extreme temperatures
          if(last!=current){# Se completo medio round trip
            # print(paste0("simulation ",ip[i,1]," iteration ",i," updating ",last," to ",current));
            rt_count[ip[i,1],c-1]=rt_count[ip[i,1],c-1]+0.5;#añadir medio round trip
            last <- current;#actualizar el current
          }  
        } 
      }
    }
  }

  ip |> filter(V4==1|V4==4) |> select(sim,V4) |> filter(sim==1)
  ip |> filter(V1==1|V1==4) |> select(sim,V1) |> filter(sim==1)
  ip |> filter(V2==1|V2==4) |> select(sim,V2) |> filter(sim==1)
  ip |> filter(V3==1|V3==4) |> select(sim,V3) |> filter(sim==1)
  ip |> filter(V3==1|V3==4) |> select(sim,V3) |> filter(sim==9)
  ip |> filter(V1==1|V1==4) |> select(sim,V1) |> filter(sim==9)
  ip|> filter(sim==9) |> pull(V1)
  ip |> filter(sim==9) |> pull(V3)
  
  return(rt_count);
}
##### Low-dimensional multimodal example #####
Rcpp::sourceCpp("functions/cpp_functions.cpp")

##### simulations for multimodal example #####
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
plot(pi.true)
plot(pi.true[-(modes+1)])

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


# Total Variation Distance
apply(output[["est_pi"]], 2,TVD,pi.est=pi.true)

# Time of first visit
output[["visits"]][modes+1]


#round trip rate
PT_RT(output[["ip"]], total_iter,iterswap,total_simulations)

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
