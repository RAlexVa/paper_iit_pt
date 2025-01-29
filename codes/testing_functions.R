#if(!require('Rcpp')){install.packages('Rcpp')}
library(Rcpp)
library(RcppArmadillo)
# setwd('..')
Rcpp::sourceCpp("functions/cpp_functions.cpp")


##### Testing IPT #####
library(Rcpp)
library(RcppArmadillo)
# setwd('..')
Rcpp::sourceCpp("functions/cpp_functions.cpp")
temperature <- c(1,0.18,0.09,.001)
bal_f <- c("sq","sq","sq","sq")
p <- 16
tot_rep <- length(temperature)

set.seed(432)
state_matrix <- matrix(rbinom(tot_rep*p,size=1,prob=0.3),ncol=tot_rep,nrow=p)
index_process <- c(3,1,4,2)
weight_matrix <- matrix(NA,nrow=p+1,ncol=tot_rep)

### Build weight matrix
for(r in 1:tot_rep){
  # print(index_process)
  currentX <- state_matrix[,r]
  logpi_current <- loglik(currentX);
  current_temp <- temperature[index_process[r]];
  chosen_bf <- bal_f[index_process[r]];
  for(n in 1:p){
    newX <-  currentX;
    newX[n] <- 1-newX[n];
    temporal <- loglik(newX)-logpi_current;
    weight_matrix[n,r] <- bal_func(temporal*current_temp, chosen_bf);
  }
  
  if(index_process[r]==4){#Last temperature
    temp_to <- temperature[1]
    index_to <- which(index_process==1)
  }else{
    temp_to <- temperature[index_process[r]+1]
    index_to <- which(index_process==(index_process[r]+1))
  }
  weight_matrix[p+1,r] <- (temp_to-current_temp)*(logpi_current - loglik(state_matrix[,index_to]))
}

##### Test function
# set.seed(43)
IPT_update(weight_matrix,state_matrix,temperature,index_process-1,bal_f)


##### #####

TVD <- function(pi,pi.est){
  return(0.5*sum(abs(pi-pi.est)))
}

X <- c(0,0,0,0,0)

IIT_update_w(X,"sq",1)

a_IIT_update(X,"sq",1,0)



### Testing 
processStrings(c("sq","min","min","min","min","sq","sq","sq","sq"))
testassignment(5, c(1,.5,.3),c("sq","sq","sq"))


#Maximo tenemos hasta dimension 30, pero es suficiente con eso
vec_to_num(c(rep(1,30),1))
vec_to_num(c(rep(0,31),1))


######## Simple example
#For p=2
set.seed(123)
p <- 8 #dimension
total_iter <- 200000
# (int p,int startsim,int endsim, int numiter,int iterswap, vec temp, const std::vector<std::string>& bal_function, bool bias_fix)
ex1 <- PT_IIT_sim(p,startsim=1, endsim=1,numiter=total_iter,iterswap=1000,temp=c(1,0.3,0.05),bal_function=c("sq","sq","sq"), bias_fix = TRUE)
ex2 <- PT_IIT_sim(p,startsim=1, endsim=1,numiter=total_iter,iterswap=1000,temp=c(1,0.3,0.05),bal_function=c("min","min","min"), bias_fix = TRUE)
est <- ex1[[1]]
est2 <- ex2[[1]]

pi1 <- est/colSums(est)
pi2 <- est2/colSums(est2)
TVD(pi1,pi2)

# Check replica swap acceptance
ip1 <- ex1[[2]]
ip2 <- ex2[[2]]

swap_rate1 <- colSums(ip1[-1, ] != ip1[-nrow(ip1), ])
swap_rate2 <- colSums(ip2[-1, ] != ip2[-nrow(ip2), ])

swap_rate1/(nrow(ip1)-1)
swap_rate2/(nrow(ip2)-1)

(constant <- (1+exp(-1))^p)

c(exp(0),exp(-1),exp(-1),exp(-2))/constant


#With no swap



loglik(c(0,0));loglik(c(1,0));loglik(c(0,1));loglik(c(1,1))

exp(c(loglik(c(0,0)),loglik(c(1,0)),loglik(c(0,1)),loglik(c(1,1))))/sum(exp(c(loglik(c(0,0)),loglik(c(1,0)),loglik(c(0,1)),loglik(c(1,1)))))



IIT_update_w(c(0,0),"sq",1)$Z
IIT_update_w(c(1,0),"sq",1)$Z
IIT_update_w(c(0,1),"sq",1)$Z
IIT_update_w(c(1,1),"sq",1)$Z


IIT_update_w(c(0,0),"min",1)$Z
IIT_update_w(c(1,0),"min",1)$Z
IIT_update_w(c(0,1),"min",1)$Z
IIT_update_w(c(1,1),"min",1)$Z


#For p=3
p <- 3 #dimension
# (int p,int startsim,int endsim, int numiter,int iterswap, vec temp, const std::vector<std::string>& bal_function, bool bias_fix)
ex1 <- PT_IIT_sim(p,startsim=1, endsim=5,numiter=50000,iterswap=100,temp=c(1,0.7,0.5),bal_function=c("sq","sq","sq"), bias_fix = TRUE)
est <- ex1[[1]]
est/colSums(est)

(constant <- (1+exp(-1))^p)
exp(0)+ 3*exp(-1) + 3*exp(-2) + exp(-3)


a <- est[,5]/sum(est[,5])

IIT_update_w(rep(0,8),"sq",1)$Z
IIT_update_w(rep(1,8),"sq",1)$Z

IIT_update_w(c(1,0,0,1,0,1,0,1),"sq",1)$Z
IIT_update_w(c(1,1,1,1,0,0,0,0),"sq",1)$Z




#probability of the mode
1/(1+exp(-1))^p

est[1,1]/sum(est[,1])

colSums(est)


#if(!require('Rcpp')){install.packages('Rcpp')}
library(Rcpp)
library(RcppArmadillo)
# setwd('..')
Rcpp::sourceCpp("functions/testing_cpp_functions.cpp")

a <- num_to_vec(16,16)
a <- num_to_vec(2^16 -1,16)
a <- num_to_vec(2^16,16)
a <- num_to_vec(0,16)
as.numeric(t(a))

as.numeric(t(num_to_vec(modes[7],16)))


sum_cube()

neigh <- matrix(c(1:5,30),nrow=2)
set.seed(123)
sample_prop(neigh)


entries_vec(0,c(1,2,3,4,5,6))
entries_vec(5,c(1,2,3,4,5,6))

##### Testing High dimensional matrix ##### 

library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp")

a <- readSparseMatrix(file.path(getwd(),"gset/G_example.txt"))



testing_loglik("gset/G_example.txt",rep(0,10))

source("functions/r_functions.R")
M <- readMatrix("gset/G_example.txt")
dim <- nrow(M)
v <- rep(0,dim)

for(i in 1:length(v)){
  tempv <- v
  tempv[i] <- 1-tempv[i]
  # print(sum(tempv))
 print(tempv%*%M%*%tempv) 
}

