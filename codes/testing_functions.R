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
# state_matrix <- matrix(rbinom(tot_rep*p,size=1,prob=0.3),ncol=tot_rep,nrow=p)
state_matrix <- matrix(0,ncol=tot_rep,nrow=p)
# index_process <- c(3,1,4,2)
index_process <- 1:4;
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
    weight_matrix[n,r] <- bal_func(temporal*current_temp, chosen_bf)-log(p);
  }
  
  if(index_process[r]==4){#Last temperature
    temp_to <- temperature[1]
    index_to <- which(index_process==1)
  }else{
    temp_to <- temperature[index_process[r]+1]
    index_to <- which(index_process==(index_process[r]+1))
  }
  temporal <- (temp_to-current_temp)*(logpi_current - loglik(state_matrix[,index_to]))
  weight_matrix[p+1,r] <-bal_func(temporal, chosen_bf)-log(tot_rep);
}

##### Test function
# set.seed(43)
ip <- index_process-1

#manual input for testing
# weight_matrix[p+1,3] <- 30
# weight_matrix[3,1] <- 40
IPT_update(weight_matrix,state_matrix,temperature,ip,bal_f)


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
M <- readMatrix("gset/ex_paper.txt")
x <- c(0,1,1,0,0)

x%*%M%*%x

y <- c(1,0,0,1,0)

y%*%M%*%y

dim <- nrow(M)
v <- rep(0,dim)

for(i in 1:length(v)){
  tempv <- v
  tempv[i] <- 1-tempv[i]
  # print(sum(tempv))
 print(tempv%*%M%*%tempv) 
}


testing_loglik("gset/ex_paper.txt",x)
print_matrix("gset/ex_paper.txt")


Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp")
# Testing matrix of size 800
print_matrix("gset/G1.txt")
v <- rep(1,800)
set.seed(345)
v2 <- rbinom(800,size=1,prob=0.4)
v3 <- rep(0:1,400)

check <- testing_loglik("gset/G1.txt",v)

######### testing updating functions
Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp")

check <- testing_loglik("gset/G1.txt",rep(0,800))

testing_a_IIT_update("gset/G1.txt",rep(0,800),"sq",1,0)



##### Testing minimum #####
library(Rcpp)
library(RcppArmadillo)
# setwd('..')
Rcpp::sourceCpp("functions/testing_cpp_functions.cpp")

x <- c(0,0,0,0,0)

find_min(x,5);x;
find_min(x,4);x;
find_min(x,2);x;
find_min(x,1);x;
find_min(x,40);x;
find_min(x,65);x;
find_min(x,-1);x;
find_min(x,100);x
##### Testing geometric distribution
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("functions/testing_cpp_functions.cpp")

print_geom(0.5)


##### Testing simualtions for high dimensional
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp")
p <- 800
temperature <- c(1,0.6,0.4)
bal_f <- c("sq","sq","sq")
set.seed(123)
#PT_a_IIT_sim(int p,int startsim,int endsim, int total_swaps,int sample_inter_swap,int burn_in, vec temp, const std::vector<std::string>& bal_function,const std::string& filename,int num_states_visited)
check3 <- PT_a_IIT_sim(p,1,2,5,100,200,temperature,bal_f,"gset/G1.txt",20)

#PT_IIT_sim(int p,int startsim,int endsim, int numiter, int iterswap,int burn_in, vec temp, const std::vector<std::string>& bal_function, bool bias_fix,const std::string& filename,int num_states_visited)
set.seed(123)
check <- PT_IIT_sim(p,1,2,500,100,200,temperature,bal_f,TRUE,"gset/G1.txt",20)
#(int p,int startsim,int endsim, int numiter,int iterswap, vec temp, const std::vector<std::string>& bal_function, bool bias_fix,const std::string& filename,int num_states_visited)
set.seed(123)
check2 <- PT_IIT_sim(p,1,2,500,100,200,temperature,bal_f,FALSE,"gset/G1.txt",20)

##### Testing burn-in period for low dimensional
rm(list=ls())
source("functions/r_functions.R")
Rcpp::sourceCpp("functions/cpp_functions.cpp")

temperature <- c(1,0.18,0.09,.001)
bal_f <- c("sq","sq","sq","sq")
p <- 16

#PT_IIT_sim(int p,int startsim,int endsim, int numiter,int iterswap,int burn_in, vec temp, const std::vector<std::string>& bal_function, bool bias_fix, int initial_state)
set.seed(123)
check <- PT_IIT_sim(p,1,2,10000,1000,5000,temperature,bal_f,TRUE,20)

#PT_IIT_sim(int p,int startsim,int endsim, int numiter,int iterswap,int burn_in, vec temp, const std::vector<std::string>& bal_function, bool bias_fix, int initial_state)
set.seed(123)
check2 <- PT_IIT_sim(p,1,2,10000,1000,5000,temperature,bal_f,FALSE,20)

set.seed(123)
#PT_a_IIT_sim(int p,int startsim,int endsim, int total_swaps,int sample_inter_swap,int burn_in, vec temp, const std::vector<std::string>& bal_function, int initial_state)
check3 <- PT_a_IIT_sim(p,1,2,10,1000,5000,temperature,bal_f,20)

set.seed(123)
#PT_a_IIT_sim(int p,int startsim,int endsim, int total_swaps,int sample_inter_swap,int burn_in, vec temp, const std::vector<std::string>& bal_function, int initial_state)
check4 <- PT_a_IIT_sim_RF(p,1,2,10000,1000,5000,temperature,bal_f,TRUE,20)





# id	algorithm	simulations	iterations	interswap	total_swap	start_state	seed	bf1	bf2	bf3	t1	t2	t3
# 4	PT_A_IIT	50	NA	1000	200	0	123	sq	sq	sq	1	0.18	0.09
rm(list=ls())
source("functions/r_functions.R")
Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp")
id_chosen <- 4
p <- 800
temperature <- c(1,0.18,0.09)
bal_f <- c("sq","sq","sq")
total_simulations <- 50
total_swap <- 200
set.seed(123)

output <- PT_a_IIT_sim_RF(p,startsim=1, endsim=2,numiter=100,iterswap=25,burn_in=100,temperature,bal_f,TRUE,"gset/G1.txt",20)


output <- PT_a_IIT_sim(p,1,50,200,1000,temperature,bal_f,"gset/G1.txt",30)

output <- PT_a_IIT_sim(p,startsim=1, endsim=total_simulations,total_swap,sample_inter_swap,temperatures,bal_f)


export <- list();
#Number of iterations needed between swaps for each replica
export[["total_iter"]] <- output[["total_iter"]]
#round trip rate (NA for IIT)
export[["round_trips"]] <- PT_RT(output[["ip"]],total_swap,total_simulations)
export[["swap_rate"]] <- output[["swap_rate"]]



export[["states"]] <- output[["states"]]
export[["loglik_visited"]] <- output[["loglik_visited"]]
export[["iter_visit"]]<- output[["iter_visit"]]
output_name <- paste0("sim_highdim_id_",id_chosen,".Rds")
saveRDS(export,file=file.path("results",output_name))



##### Testing binary vector #####
rm(list=ls())
source("functions/r_functions.R")
Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp")

a <- createBinaryVector(c(1,500,800,350,850,-1),800)



##### Checking swap rate #####
rm(list=ls())
source("functions/r_functions.R")
Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp")
id_chosen <- 4
p <- 800
temperature <- c(1,0.18,0.09)
bal_f <- c("sq","sq","sq")
total_simulations <- 50
total_swap <- 200
set.seed(123)

output <- PT_IIT_sim(p,startsim=1, endsim=2,numiter=100,iterswap=25,burn_in=100,temperature,bal_f,TRUE,"gset/G1.txt",20)


#### Testing initializing the high dimensional problem
rm(list=ls())
Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp")
createBinaryVector(c(1,4,5,10)-1,10)

initializeMatrix(c(1,4,5,10)-1,10,5)



#### Testing random likelihoods
rm(list=ls())
Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp")
check <- hist_lik("gset/G1.txt",500000,0.3)
hist(check)
max(check)
min(check)
check <- hist_lik("gset/G1.txt",500000,0.7)
hist(check)
max(check)
min(check)



#### Testing likelihoods of multimodal highdimensional
rm(list=ls())
Rcpp::sourceCpp("functions/cpp_func_multihigh.cpp")
vec1 <- c(rep(0,400),rep(1,400))
vec2 <- c(rep(1,400),rep(0,400))
vec3 <- rep(0:1,400)

eval_loglik("gset/G1.txt",vec1)
eval_loglik("gset/G1.txt",vec2)
eval_loglik("gset/G1.txt",vec3)



#### Testing bounded balancing function
rm(list=ls())
Rcpp::sourceCpp("functions/cpp_functions.cpp")
bound_sq(1,2)

#Using bound 2 we get 
library(ggplot2)
exp(bound_sq(log(0.1),log(2)))
exp(bound_sq(log(0.2),log(2)))
exp(bound_sq(log(0.25),log(2)))
exp(bound_sq(log(0.5),log(2)))
exp(bound_sq(log(3),log(2)))
exp(bound_sq(log(54),log(2)))
exp(bound_sq(log(55),log(2)))
exp(bound_sq(log(155),log(2)))
exp(bound_sq(log(12355),log(2)))
y <- c()
x <- c(seq(0.01,0.25,by=0.01),seq(0.25,5,by=0.1))
x <- c(seq(0.01,0.25,by=0.001))
for(i in 1:length(x)){
  y[i] <- exp(bound_sq(log(x[i]),log(2)))
}

as.data.frame(x=x,y=y) |> ggplot(aes(x=x,y=y))+geom_point()



######
##### Testing adaptivve PT (decreasing of bounding constant)
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("functions/cpp_functions.cpp")
p <- 16
temperature <- c(1,0.155,0.1,0.05)
bal_f <- c("sq","sq","sq","sq")
set.seed(123)
#PT_a_IIT_sim(p,1, total_simulations,total_swap,sample_inter_swap,burnin_iter,temperatures,bal_f,start_state,apply_reduction,reduc_constant,reduc_model)
check <- PT_a_IIT_sim(p,1,1,total_swap=100,sample_inter_swap=10000,200,temperature,bal_f,0,TRUE,.1,"iterations")


#Testing optimum of G1

op1 <- c(0,1,1,1,0,0,0,0,0,1,1,0,0,1,0,1,1,0,0,0,1,1,0,0,0,0,1,0,1,0,1,1,1,0,0,1,0,0,1,1,0,1,1,1,1,1,1,1,0,1,0,1,0,1,0,1,0,1,0,0,1,1,0,0,0,0,0,0,1,1,0,0,1,0,1,1,1,1,0,0,0,1,1,1,0,1,1,1,0,0,0,1,0,1,1,1,1,1,0,1,1,0,1,1,0,1,0,0,1,1,0,1,1,1,0,0,0,0,1,1,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,1,0,1,1,0,0,1,1,1,1,1,0,1,1,0,1,0,1,0,1,0,1,1,0,1,0,1,1,0,0,1,0,1,0,1,1,0,1,1,1,0,0,1,1,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0,1,0,0,0,1,1,1,0,0,0,1,0,0,1,0,1,1,1,1,0,0,0,1,0,0,0,1,1,1,1,0,1,0,1,0,1,1,1,1,0,0,1,1,0,0,0,0,0,0,0,1,1,1,1,1,0,1,1,1,1,1,1,1,0,0,0,1,0,0,1,0,0,0,1,0,0,0,1,1,1,1,1,0,1,0,1,1,1,1,1,0,1,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,1,1,1,0,0,1,1,0,1,0,0,0,1,1,0,0,1,1,0,1,1,1,0,0,0,1,1,1,1,0,0,1,1,1,1,0,1,0,0,1,1,0,0,0,1,0,1,1,0,0,1,0,1,1,0,1,0,0,1,0,1,0,1,1,0,1,0,0,1,0,0,1,0,1,0,1,1,1,0,1,0,0,1,0,1,0,1,1,0,1,0,1,0,1,1,0,0,1,1,0,1,0,0,0,1,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,1,0,1,0,1,1,1,0,1,0,0,0,1,0,1,0,0,1,1,0,1,0,1,1,0,1,1,0,0,1,1,1,1,1,0,0,1,0,0,0,1,1,0,0,1,1,0,0,0,1,1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0,1,1,1,0,0,0,1,0,0,1,1,1,1,0,1,1,1,1,1,1,0,0,0,1,1,1,1,1,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1,0,1,0,0,1,1,0,1,0,1,1,1,1,0,0,0,0,0,1,0,0,0,1,0,0,1,1,0,0,0,0,1,0,1,1,1,0,1,0,1,1,0,1,0,1,1,0,0,1,1,1,0,0,0,1,1,0,0,1,1,1,1,0,0,0,1,0,1,1,0,0,1,1,0,1,1,0,1,1,0,1,1,1,0,0,0,0,0,0,0,1,1,0,1,1,1,0,1,0,0,0,0,1,0,0,1,0,0,1,1,1,0,0,0,0,1,1,1,1,0,1,1,0,1,1,1,0,1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,0,1,1,1,0,1,0,1,1,1,0,1,1,1,1,1,0,1,1,0,0,1,0,0,1,0,1,0,1,0,1,0,1,0,1,0,0,0,0,1,1,1,1,0,0,1,1,0,0,0,0,0,1,0,0,1,1,1,0,1,0,0,0,1,0,0,0,0,0,1,1,0,0,1,0,0,0,1,0,0,1,0,0,1,1,0,0,1,1,0,1,1,1,1,0,1,0,0,1,1)
op2 <- c(1,0,0,0,1,1,1,1,1,0,0,1,1,0,1,0,0,1,1,1,0,0,1,1,1,1,0,1,0,1,0,0,0,1,1,0,1,1,0,0,1,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,1,0,0,1,1,1,1,1,1,0,0,1,1,0,1,1,0,0,0,1,1,1,0,0,0,1,0,0,0,1,1,1,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,1,1,0,0,1,0,0,0,1,1,1,1,0,0,1,1,1,1,1,0,0,0,1,0,1,1,1,1,1,0,1,0,0,1,1,0,0,0,0,0,1,0,0,1,0,1,0,1,0,1,0,0,1,0,1,0,0,1,1,0,1,0,1,0,0,1,0,0,0,1,1,0,0,0,1,1,1,1,0,1,0,1,1,1,1,1,0,1,0,1,1,1,0,0,0,1,1,1,0,1,1,0,1,0,0,0,0,1,1,1,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,0,1,1,1,0,1,1,1,1,1,1,1,0,0,0,0,0,1,0,0,0,1,0,0,0,1,1,1,0,1,1,0,1,1,1,0,1,1,1,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0,1,1,1,0,1,1,1,0,1,0,1,1,1,1,1,1,0,0,0,1,1,0,0,1,0,1,1,1,0,0,1,1,0,0,1,0,0,0,1,1,1,0,0,0,0,1,1,0,0,0,0,1,0,1,1,0,0,1,1,1,0,1,0,0,1,1,0,1,0,0,1,1,1,1,0,1,0,1,0,0,1,0,1,1,0,1,1,0,1,0,1,0,0,0,1,0,1,1,0,1,0,1,0,0,1,0,1,0,1,0,0,1,1,0,0,1,0,1,1,1,0,1,1,1,0,0,0,1,1,1,1,1,1,0,0,0,0,1,0,1,0,0,0,1,0,1,1,1,0,1,0,1,1,0,0,1,0,1,0,0,1,0,0,1,1,0,0,0,0,0,1,1,0,1,1,1,0,0,1,1,0,1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,1,1,1,0,0,1,1,1,0,0,1,1,1,1,0,1,1,0,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,1,1,0,1,1,1,1,1,0,1,0,1,0,1,1,0,0,1,0,1,0,0,0,0,1,1,0,1,1,0,1,1,1,0,1,1,0,0,1,1,1,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,1,1,0,0,0,1,1,1,0,0,1,1,0,0,0,0,1,1,1,0,1,0,0,1,1,0,0,1,0,0,1,0,0,1,0,0,0,1,1,1,1,1,1,1,0,0,1,0,0,0,1,0,1,1,1,1,0,1,1,0,1,1,0,0,0,1,1,1,1,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,0,0,1,0,1,1,1,1,1,1,0,1,1,1,1,1,0,0,0,1,1,1,0,1,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,1,0,0,1,1,0,1,1,0,1,1,1,0,1,0,1,0,1,0,1,1,1,1,0,0,0,0,1,1,0,0,1,1,1,1,1,0,1,1,0,0,0,1,0,1,1,1,0,1,1,1,1,1,0,0,1,1,0,1,1,1,0,1,1,0,1,0,0,0,1,1,0,0,1,0,0,0,0,1,0,1,1,0,0)
Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp")
eval_loglik("gset/G1.txt",op1)
eval_lik("gset/G1.txt",op1)
eval_loglik("gset/G1.txt",op2)
eval_lik("gset/G1.txt",op2)

s1 <- testing_lik("gset/G1.txt",op1)
s2 <- testing_lik("gset/G1.txt",op2)



####### Testing how the constant changes
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("functions/testing_cpp_functions.cpp")

n <- 1815
initial_bound <- 0.375
# print_log_bound(int iterations, double initial_bound, double prob_to_dec, double temperature, double decreasing_constant)
print_log_bound(n,initial_bound,1,.05,1)


n <- 18
initial_bound <- 0.05
# print_log_bound(int iterations, double initial_bound, double prob_to_dec, double temperature, double decreasing_constant)
print_log_bound(n,initial_bound,1,.05,1)


n <- 1815
initial_bound <- 7.5
# print_log_bound(int iterations, double initial_bound, double prob_to_dec, double temperature, double decreasing_constant)
print_log_bound(n,initial_bound,1,1,1)

n <- 1000
initial_bound <- 0.375
# print_log_bound(int iterations, double initial_bound, double prob_to_dec, double temperature, double decreasing_constant)
print_log_bound(n,initial_bound,1,.05,.1)


### Testing high values of log bound
Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp")

v1 <- rep(0,800)
v1[c(102,197,227,437,445)-1] <- 1

c_l <- eval_lik("gset/G1.txt",v1)
check <- testing_lik("gset/G1.txt",v1)
max(check)
min(check)

c_l-min(check)
max(check)-c_l


c_l <- eval_loglik("gset/G1.txt",v1)
check <- testing_loglik("gset/G1.txt",v1)
max(check)
min(check)

c_l-min(check)
max(check)-c_l



#### testing high dimensional multimodal problem
Rcpp::sourceCpp("functions/cpp_functions_highdim_2.cpp")

M <- matrix(c(1,1,1,1,0,1,0,0),ncol=2,byrow=F)
v <- c(1,0,1,1)
loglik(v,M)

mod1 <- rep(0:1,400)
mod2 <- rep(1:0,400)
M <- cbind(mod1,mod2)


mod1 <- c(rep(1,150),rep(0,150),rep(1,200),rep(0,300))
mod2 <- c(rep(0,150),rep(1,150),rep(1,200),rep(0,300))
M <- cbind(mod1,mod2)

loglik(mod1,M)
loglik(mod2,M)

loglik(rep(0,800),M)

loglik(rbinom(800,1,0.3),M)
lvec <- c()
for(i in 1:length(mod1)){
  temp_vec <- mod1
  temp_vec[i] <- 1-temp_vec[i]
  lvec[i] <- loglik(temp_vec,M)
}
summary(lvec)

lvec <- c()
for(i in 1:length(mod2)){
  temp_vec <- mod2
  temp_vec[i] <- 1-temp_vec[i]
  lvec[i] <- loglik(temp_vec,M)
}
summary(lvec)


#### testing bimodal high dim example
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("functions/cpp_functions_highdim_2.cpp")

p <- 200
Q_matrix <- matrix(0,nrow=p,ncol=2)
for(i in 1:p){
  if(i%%2==1){Q_matrix[i,2]=1}
  if(i%%2==0){Q_matrix[i,1]=1}
}

v0 <- rep(0,p)
v1 <- Q_matrix[,1]
v2 <- Q_matrix[,2]

loglik(data$states[,1,1],Q_matrix)
loglik(data$states[,4,1],Q_matrix)

loglik(v0,Q_matrix)
loglik(v1,Q_matrix)
loglik(v2,Q_matrix)
exp(-40)
log(1+exp(-40))
# loglik(c(rep(0,400),rep(1,400)),Q_matrix)

vr <- rbinom(p,1,0.5)
loglik(vr,Q_matrix)
loglik(rbinom(p,1,0.5),Q_matrix)

check <- a_IIT_update(v0,Q_matrix,"sq",1,0)

check_res <- PT_a_IIT_sim(p,1,1, 5,1,0, c(1,0.3),c("sq","sq"),"nothing",5,-1)

create_vector <- function(coord,dim){
  vec <- rep(0,dim)
  vec[coord] <- 1
  return(vec)
}

v_check <- v0
v_check[170] <- 1

loglik(v_check,Q_matrix)

store_lik <- c()
for(i in 1:length(v1)){
  rep_v <- v1
  rep_v[i] <- 1-rep_v[i]
  store_lik[i] <- loglik(rep_v,Q_matrix)
}
summary(store_lik)


neigh_1 <- v1
neigh_1[30] <- 1-neigh_1[30]
loglik(neigh_1,Q_matrix)
store_lik <- c()
for(i in 1:length(neigh_1)){
  rep_v <- neigh_1
  rep_v[i] <- 1-rep_v[i]
  store_lik[i] <- loglik(rep_v,Q_matrix)
}
summary(store_lik)
