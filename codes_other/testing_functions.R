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
Rcpp::sourceCpp("functions_other/testing_cpp_functions.cpp")

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

### Tracing a path to the modes
create_vector <- function(coord,dim){
  vec <- rep(0,dim)
  vec[coord] <- 1
  return(vec)
}
path <- c()
for(i in 1:p){
  if(i%%2==1){
    path[(i+1)/2] <- loglik(create_vector(seq(1,i,by=2),p),Q_matrix) 
  }
}



exp(-40)
log(1+exp(-40))
# loglik(c(rep(0,400),rep(1,400)),Q_matrix)

vr <- rbinom(p,1,0.5)
loglik(vr,Q_matrix)
loglik(rbinom(p,1,0.5),Q_matrix)

check <- a_IIT_update(v0,Q_matrix,"sq",1,0)

check_res <- PT_a_IIT_sim(p,1,1, 5,1,0, c(1,0.3),c("sq","sq"),"nothing",5,-1)



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


### Testing time of first visit
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("functions/cpp_functions.cpp")
p <- 16

#PT_IIT_sim(int p,int startsim,int endsim, int numiter,int iterswap,int burn_in, vec temp, const std::vector<std::string>& bal_function, bool bias_fix, int initial_state)
check_res <- PT_IIT_sim(p,1,5, 5000,100,1000, c(1,0.3),c("sq","sq"),TRUE,654)


##### Checking how to visit at least 1 mode in the high dimensional bimodal example #####
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("functions/cpp_functions_highdim_2.cpp")
p <- 200
numiter <- 50000
burn_in <- 10000


Q_matrix <- matrix(0,nrow=p,ncol=2)
for(i in 1:p){
  if(i%%2==1){Q_matrix[i,2]=1}
  if(i%%2==0){Q_matrix[i,1]=1}
}

v0 <- rep(0,p)
v1 <- Q_matrix[,1]
v2 <- Q_matrix[,2]

loglik(v0,Q_matrix)
loglik(v1,Q_matrix)
loglik(v1,Q_matrix)
loglik(rbinom(p,1,0.4),Q_matrix)

set.seed(456)
temperatures <- c(1.5,1.2,1,0.85,0.75)
bal_f <- rep("sq",length(temperatures))

#This one the 3 first replicas find both modes, the first 2 in the burn-in and temperature 1 in iteration 35k
set.seed(456)
check_res <- PT_IIT_sim(p,1,1,numiter, 10,burn_in, temperatures,bal_f, TRUE,"anything",3,-1)


#This one the 3 first replicas find both modes but it takes like... 75k iterations
set.seed(456)
check_res <- PT_IIT_sim(p,1,1,numiter, 50,burn_in, temperatures,bal_f, TRUE,"anything",3,-1)

#PT_a_IIT_sim(int p,int startsim,int endsim, int total_swaps,int sample_inter_swap,int burn_in, vec temp, const std::vector<std::string>& bal_function,const std::string& filename,int num_states_visited,const std::vector<int>& starting_coord)
p <- 200
burn_in <- 30000
total_swaps <- 250
inter_swap <- 1000
temperatures <- c(1.4,1.25,1.13,1,0.9,0.82)
bal_f <- rep("sq",length(temperatures))
set.seed(456)
check_ad <- PT_a_IIT_sim(p,1,1,total_swaps,inter_swap,burn_in, temperatures, bal_f,"anything",3,-1)
check_ad$swap_rate
check_ad$modes_visit
check_ad$time_taken

##### Tracing a path to the mode of the Gset high dim example #####

Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp")
p <- 800
op1 <- c(0,1,1,1,0,0,0,0,0,1,1,0,0,1,0,1,1,0,0,0,1,1,0,0,0,0,1,0,1,0,1,1,1,0,0,1,0,0,1,1,0,1,1,1,1,1,1,1,0,1,0,1,0,1,0,1,0,1,0,0,1,1,0,0,0,0,0,0,1,1,0,0,1,0,1,1,1,1,0,0,0,1,1,1,0,1,1,1,0,0,0,1,0,1,1,1,1,1,0,1,1,0,1,1,0,1,0,0,1,1,0,1,1,1,0,0,0,0,1,1,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,1,0,1,1,0,0,1,1,1,1,1,0,1,1,0,1,0,1,0,1,0,1,1,0,1,0,1,1,0,0,1,0,1,0,1,1,0,1,1,1,0,0,1,1,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0,1,0,0,0,1,1,1,0,0,0,1,0,0,1,0,1,1,1,1,0,0,0,1,0,0,0,1,1,1,1,0,1,0,1,0,1,1,1,1,0,0,1,1,0,0,0,0,0,0,0,1,1,1,1,1,0,1,1,1,1,1,1,1,0,0,0,1,0,0,1,0,0,0,1,0,0,0,1,1,1,1,1,0,1,0,1,1,1,1,1,0,1,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,1,1,1,0,0,1,1,0,1,0,0,0,1,1,0,0,1,1,0,1,1,1,0,0,0,1,1,1,1,0,0,1,1,1,1,0,1,0,0,1,1,0,0,0,1,0,1,1,0,0,1,0,1,1,0,1,0,0,1,0,1,0,1,1,0,1,0,0,1,0,0,1,0,1,0,1,1,1,0,1,0,0,1,0,1,0,1,1,0,1,0,1,0,1,1,0,0,1,1,0,1,0,0,0,1,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,1,0,1,0,1,1,1,0,1,0,0,0,1,0,1,0,0,1,1,0,1,0,1,1,0,1,1,0,0,1,1,1,1,1,0,0,1,0,0,0,1,1,0,0,1,1,0,0,0,1,1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0,1,1,1,0,0,0,1,0,0,1,1,1,1,0,1,1,1,1,1,1,0,0,0,1,1,1,1,1,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1,0,1,0,0,1,1,0,1,0,1,1,1,1,0,0,0,0,0,1,0,0,0,1,0,0,1,1,0,0,0,0,1,0,1,1,1,0,1,0,1,1,0,1,0,1,1,0,0,1,1,1,0,0,0,1,1,0,0,1,1,1,1,0,0,0,1,0,1,1,0,0,1,1,0,1,1,0,1,1,0,1,1,1,0,0,0,0,0,0,0,1,1,0,1,1,1,0,1,0,0,0,0,1,0,0,1,0,0,1,1,1,0,0,0,0,1,1,1,1,0,1,1,0,1,1,1,0,1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,0,1,1,1,0,1,0,1,1,1,0,1,1,1,1,1,0,1,1,0,0,1,0,0,1,0,1,0,1,0,1,0,1,0,1,0,0,0,0,1,1,1,1,0,0,1,1,0,0,0,0,0,1,0,0,1,1,1,0,1,0,0,0,1,0,0,0,0,0,1,1,0,0,1,0,0,0,1,0,0,1,0,0,1,1,0,0,1,1,0,1,1,1,1,0,1,0,0,1,1)
op2 <- c(1,0,0,0,1,1,1,1,1,0,0,1,1,0,1,0,0,1,1,1,0,0,1,1,1,1,0,1,0,1,0,0,0,1,1,0,1,1,0,0,1,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,1,0,0,1,1,1,1,1,1,0,0,1,1,0,1,1,0,0,0,1,1,1,0,0,0,1,0,0,0,1,1,1,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,1,1,0,0,1,0,0,0,1,1,1,1,0,0,1,1,1,1,1,0,0,0,1,0,1,1,1,1,1,0,1,0,0,1,1,0,0,0,0,0,1,0,0,1,0,1,0,1,0,1,0,0,1,0,1,0,0,1,1,0,1,0,1,0,0,1,0,0,0,1,1,0,0,0,1,1,1,1,0,1,0,1,1,1,1,1,0,1,0,1,1,1,0,0,0,1,1,1,0,1,1,0,1,0,0,0,0,1,1,1,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,0,1,1,1,0,1,1,1,1,1,1,1,0,0,0,0,0,1,0,0,0,1,0,0,0,1,1,1,0,1,1,0,1,1,1,0,1,1,1,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0,1,1,1,0,1,1,1,0,1,0,1,1,1,1,1,1,0,0,0,1,1,0,0,1,0,1,1,1,0,0,1,1,0,0,1,0,0,0,1,1,1,0,0,0,0,1,1,0,0,0,0,1,0,1,1,0,0,1,1,1,0,1,0,0,1,1,0,1,0,0,1,1,1,1,0,1,0,1,0,0,1,0,1,1,0,1,1,0,1,0,1,0,0,0,1,0,1,1,0,1,0,1,0,0,1,0,1,0,1,0,0,1,1,0,0,1,0,1,1,1,0,1,1,1,0,0,0,1,1,1,1,1,1,0,0,0,0,1,0,1,0,0,0,1,0,1,1,1,0,1,0,1,1,0,0,1,0,1,0,0,1,0,0,1,1,0,0,0,0,0,1,1,0,1,1,1,0,0,1,1,0,1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,1,1,1,0,0,1,1,1,0,0,1,1,1,1,0,1,1,0,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,1,1,0,1,1,1,1,1,0,1,0,1,0,1,1,0,0,1,0,1,0,0,0,0,1,1,0,1,1,0,1,1,1,0,1,1,0,0,1,1,1,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,1,1,0,0,0,1,1,1,0,0,1,1,0,0,0,0,1,1,1,0,1,0,0,1,1,0,0,1,0,0,1,0,0,1,0,0,0,1,1,1,1,1,1,1,0,0,1,0,0,0,1,0,1,1,1,1,0,1,1,0,1,1,0,0,0,1,1,1,1,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,0,0,1,0,1,1,1,1,1,1,0,1,1,1,1,1,0,0,0,1,1,1,0,1,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,1,0,0,1,1,0,1,1,0,1,1,1,0,1,0,1,0,1,0,1,1,1,1,0,0,0,0,1,1,0,0,1,1,1,1,1,0,1,1,0,0,0,1,0,1,1,1,0,1,1,1,1,1,0,0,1,1,0,1,1,1,0,1,1,0,1,0,0,0,1,1,0,0,1,0,0,0,0,1,0,1,1,0,0)
eval_lik("gset/G1.txt",op1)
eval_lik("gset/G1.txt",1-op1)
eval_lik("gset/G1.txt",op2)


path1 <- lik_path("gset/G1.txt", rep(0,p), op1)
plot(1:800,path1)
min(path1[-800]-path1[-1])
path2 <- lik_path("gset/G1.txt", rep(0,p), op2)
plot(1:800,path2)
min(path2[-800]-path2[-1])

##### Tracing a path between modes of the Gset high dim example #####
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp")
p <- 800
op1 <- c(0,1,1,1,0,0,0,0,0,1,1,0,0,1,0,1,1,0,0,0,1,1,0,0,0,0,1,0,1,0,1,1,1,0,0,1,0,0,1,1,0,1,1,1,1,1,1,1,0,1,0,1,0,1,0,1,0,1,0,0,1,1,0,0,0,0,0,0,1,1,0,0,1,0,1,1,1,1,0,0,0,1,1,1,0,1,1,1,0,0,0,1,0,1,1,1,1,1,0,1,1,0,1,1,0,1,0,0,1,1,0,1,1,1,0,0,0,0,1,1,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,1,0,1,1,0,0,1,1,1,1,1,0,1,1,0,1,0,1,0,1,0,1,1,0,1,0,1,1,0,0,1,0,1,0,1,1,0,1,1,1,0,0,1,1,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0,1,0,0,0,1,1,1,0,0,0,1,0,0,1,0,1,1,1,1,0,0,0,1,0,0,0,1,1,1,1,0,1,0,1,0,1,1,1,1,0,0,1,1,0,0,0,0,0,0,0,1,1,1,1,1,0,1,1,1,1,1,1,1,0,0,0,1,0,0,1,0,0,0,1,0,0,0,1,1,1,1,1,0,1,0,1,1,1,1,1,0,1,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,1,1,1,0,0,1,1,0,1,0,0,0,1,1,0,0,1,1,0,1,1,1,0,0,0,1,1,1,1,0,0,1,1,1,1,0,1,0,0,1,1,0,0,0,1,0,1,1,0,0,1,0,1,1,0,1,0,0,1,0,1,0,1,1,0,1,0,0,1,0,0,1,0,1,0,1,1,1,0,1,0,0,1,0,1,0,1,1,0,1,0,1,0,1,1,0,0,1,1,0,1,0,0,0,1,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,1,0,1,0,1,1,1,0,1,0,0,0,1,0,1,0,0,1,1,0,1,0,1,1,0,1,1,0,0,1,1,1,1,1,0,0,1,0,0,0,1,1,0,0,1,1,0,0,0,1,1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0,1,1,1,0,0,0,1,0,0,1,1,1,1,0,1,1,1,1,1,1,0,0,0,1,1,1,1,1,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1,0,1,0,0,1,1,0,1,0,1,1,1,1,0,0,0,0,0,1,0,0,0,1,0,0,1,1,0,0,0,0,1,0,1,1,1,0,1,0,1,1,0,1,0,1,1,0,0,1,1,1,0,0,0,1,1,0,0,1,1,1,1,0,0,0,1,0,1,1,0,0,1,1,0,1,1,0,1,1,0,1,1,1,0,0,0,0,0,0,0,1,1,0,1,1,1,0,1,0,0,0,0,1,0,0,1,0,0,1,1,1,0,0,0,0,1,1,1,1,0,1,1,0,1,1,1,0,1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,0,1,1,1,0,1,0,1,1,1,0,1,1,1,1,1,0,1,1,0,0,1,0,0,1,0,1,0,1,0,1,0,1,0,1,0,0,0,0,1,1,1,1,0,0,1,1,0,0,0,0,0,1,0,0,1,1,1,0,1,0,0,0,1,0,0,0,0,0,1,1,0,0,1,0,0,0,1,0,0,1,0,0,1,1,0,0,1,1,0,1,1,1,1,0,1,0,0,1,1)
op2 <- c(1,0,0,0,1,1,1,1,1,0,0,1,1,0,1,0,0,1,1,1,0,0,1,1,1,1,0,1,0,1,0,0,0,1,1,0,1,1,0,0,1,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,1,0,0,1,1,1,1,1,1,0,0,1,1,0,1,1,0,0,0,1,1,1,0,0,0,1,0,0,0,1,1,1,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,1,1,0,0,1,0,0,0,1,1,1,1,0,0,1,1,1,1,1,0,0,0,1,0,1,1,1,1,1,0,1,0,0,1,1,0,0,0,0,0,1,0,0,1,0,1,0,1,0,1,0,0,1,0,1,0,0,1,1,0,1,0,1,0,0,1,0,0,0,1,1,0,0,0,1,1,1,1,0,1,0,1,1,1,1,1,0,1,0,1,1,1,0,0,0,1,1,1,0,1,1,0,1,0,0,0,0,1,1,1,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,0,1,1,1,0,1,1,1,1,1,1,1,0,0,0,0,0,1,0,0,0,1,0,0,0,1,1,1,0,1,1,0,1,1,1,0,1,1,1,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0,1,1,1,0,1,1,1,0,1,0,1,1,1,1,1,1,0,0,0,1,1,0,0,1,0,1,1,1,0,0,1,1,0,0,1,0,0,0,1,1,1,0,0,0,0,1,1,0,0,0,0,1,0,1,1,0,0,1,1,1,0,1,0,0,1,1,0,1,0,0,1,1,1,1,0,1,0,1,0,0,1,0,1,1,0,1,1,0,1,0,1,0,0,0,1,0,1,1,0,1,0,1,0,0,1,0,1,0,1,0,0,1,1,0,0,1,0,1,1,1,0,1,1,1,0,0,0,1,1,1,1,1,1,0,0,0,0,1,0,1,0,0,0,1,0,1,1,1,0,1,0,1,1,0,0,1,0,1,0,0,1,0,0,1,1,0,0,0,0,0,1,1,0,1,1,1,0,0,1,1,0,1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,1,1,1,0,0,1,1,1,0,0,1,1,1,1,0,1,1,0,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,1,1,0,1,1,1,1,1,0,1,0,1,0,1,1,0,0,1,0,1,0,0,0,0,1,1,0,1,1,0,1,1,1,0,1,1,0,0,1,1,1,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,1,1,0,0,0,1,1,1,0,0,1,1,0,0,0,0,1,1,1,0,1,0,0,1,1,0,0,1,0,0,1,0,0,1,0,0,0,1,1,1,1,1,1,1,0,0,1,0,0,0,1,0,1,1,1,1,0,1,1,0,1,1,0,0,0,1,1,1,1,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,0,0,1,0,1,1,1,1,1,1,0,1,1,1,1,1,0,0,0,1,1,1,0,1,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,1,0,0,1,1,0,1,1,0,1,1,1,0,1,0,1,0,1,0,1,1,1,1,0,0,0,0,1,1,0,0,1,1,1,1,1,0,1,1,0,0,0,1,0,1,1,1,0,1,1,1,1,1,0,0,1,1,0,1,1,1,0,1,1,0,1,0,0,0,1,1,0,0,1,0,0,0,0,1,0,1,1,0,0)
(l1 <- eval_lik("gset/G1.txt",op1))
eval_lik("gset/G1.txt",1-op1)
(l2 <- eval_lik("gset/G1.txt",op2))

(lo <- eval_lik("gset/G1.txt",c(1,rep(0,p-1))))

path <- lik_path("gset/G1.txt",op1,op2)

### Definition of minimum temperature
#### Considering this we can take temperature 0.005 as a starting point for this problem
min_temperature <- 0.0005
min(path)^min_temperature
l1^min_temperature
l2^min_temperature
lo^min_temperature

prop_temp <- 0.0005
total_temps <- 20
temp_ladder <- c(1,0.95,0.90,0.8,0.7,0.5,0.4,0.3,0.2,0.1,0.09,0.05,0.01,0.0005)
total_temps <- length(temp_ladder)
temp_ladder
# temp_ladder <- c(0.0005,0.0006)
data_plot <- tibble(x=numeric(0),y=numeric(0),logy=numeric(0),temp=numeric(0))
for(i in 1:total_temps){
  x <- 1:p
  logy <- log(path)*temp_ladder[i]
  y<- path^temp_ladder[i]
  temp <- temp_ladder[i]
  data_plot <- rbind(data_plot,cbind(x,y,logy,temp))
}
colnames(data_plot) <- c("x","y","logy","temp")
data_plot |> ggplot(aes(x=x,y=y,col=temp))+geom_point(size=0.5)

data_plot |> ggplot(aes(x=x,y=logy,col=temp))+geom_point(size=0.5)



##### Checking how to visit at least 1 mode in the high dimensional bimodal example #####
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
Rcpp::sourceCpp("functions/cpp_functions_highdim_2.cpp")
p <- 200

v1 <- rep(c(0,1),p/2)
v2 <- rep(c(1,0),p/2)
Q_matrix <- cbind(v1,v2)

v0 <- rep(0,p)
loglik(v0,Q_matrix)
loglik(v1,Q_matrix)
loglik(v2,Q_matrix)
loglik(rbinom(p,1,0.4),Q_matrix)

vector1 <- rep(0,p)
vector2 <- rep(0,p)
vector3 <- v1

lik_p1 <- c()
lik_p2 <- c()
lik_modes <- c()
### Tracing a path from 0 to modes and from modes to modes
for(i in 1:p){
  vector1[i] <- v1[i]
  vector2[i] <- v2[i]
  vector3[i] <- v2[i]
  lik_p1[i] <- loglik(vector1,Q_matrix)
  lik_p2[i] <- loglik(vector2,Q_matrix)
  lik_modes[i] <- loglik(vector3,Q_matrix)
}

### Definition of minimum temperature
#### Considering this we can take temperature 0.0001 as a starting point for this problem
min_temperature <- 0.0001
exp(min(lik_modes)*min_temperature)
exp(loglik(v1,Q_matrix)*min_temperature)
exp(loglik(v2,Q_matrix)*min_temperature)
exp(loglik(v0,Q_matrix)*min_temperature)



plot(1:p,lik_p1)
plot(1:p,lik_p2)

dif_modes <- lik_modes[-p]-lik_modes[-1]

# prop_temp <- 0.0006
prop_temp <- 0.0001

tibble(x=1:p,y1=lik_modes, y2=lik_modes*prop_temp) |> 
  ggplot(aes(x=x,y=y1))+
  geom_line(col="blue")+
  geom_line(aes(y=y2), col="red")+
  labs(y='log-likelihood',x="step")

tibble(x=1:p,y1=exp(lik_modes), y2=exp(lik_modes*prop_temp)) |> 
  ggplot(aes(x=x,y=y1))+
  geom_line(col="blue")+
  geom_line(aes(y=y2), col="red")+
  labs(y='likelihood',x="step")

total_temps <- 20
temp_ladder <- (seq(log(prop_temp),0,length.out=total_temps))
temp_ladder
temp_ladder <- exp(temp_ladder)
temp_ladder
data_plot <- tibble(x=numeric(0),y=numeric(0),logy=numeric(0),temp=numeric(0))
for(i in 1:total_temps){
  x <- 1:p
  logy <- lik_modes*temp_ladder[i]
  y<- exp(lik_modes*temp_ladder[i])
  temp <- temp_ladder[i]
  data_plot <- rbind(data_plot,cbind(x,y,logy,temp))
}

data_plot |> ggplot(aes(x=x,y=y,col=temp))+geom_point(size=0.5)

data_plot |> ggplot(aes(x=x,y=logy,col=temp))+geom_point(size=0.5)




temperatures <- c(1.15,1.1,1.05,1,0.95,0.9,0.85)
bal_f <- rep("sq",length(temperatures))

#This one the 3 first replicas find both modes, the first 2 in the burn-in and temperature 1 in iteration 35k
set.seed(456)
numiter <- 100000
burn_in <- 1000
bidim <- PT_IIT_sim(p,1,1,numiter, 100,burn_in, temperatures,bal_f, TRUE,"anything",3,-1)






########## Tracing a path to the modes in the bimodal high dim example #####
create_vector <- function(coord,dim){
  vec <- rep(0,dim)
  vec[coord] <- 1
  return(vec)
}
path <- c()

for(i in 1:p){
  if(i%%2==1){
    path[(i+1)/2] <- loglik(create_vector(seq(1,i,by=2),p),Q_matrix) 
  }
}

loglik(vec,Q_matrix)
store_lik <- c()
for(i in 1:length(vec)){
  rep_v <- vec
  rep_v[i] <- 1-rep_v[i]
  store_lik[i] <- loglik(rep_v,Q_matrix)
}
summary(store_lik)


##### Testing implementation of bound_sq #####
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp")
p <- 800
tot_swaps <- 10
inter_swap <- 500
burn_in <- 2000
temperatures <- c(1,0.9)
bal_f <- rep("sq",length(temperatures))

set.seed(89)
(sim_c <- PT_a_IIT_sim(p,1,1,tot_swaps,inter_swap,burn_in,temperatures, bal_f,"gset/G1.txt",num_states_visited=5,-1, decreasing_constant=.000000001,"iterations"))
(sim_c <- PT_a_IIT_sim(p,1,1,tot_swaps,inter_swap,burn_in,temperatures, bal_f,"gset/G1.txt",num_states_visited=5,-1, decreasing_constant=.00005,"never"))

list_bk <- list()

for(e in names(sim_c)){
  list_bk[[e]] <- sim_c[[e]]
}



#### Testing report of TVD #####
##### Testing IPT #####
library(Rcpp)
library(RcppArmadillo)
# setwd('..')
Rcpp::sourceCpp("functions/cpp_functions.cpp")
p <- 16
total_iter <- 30000
iter_swap <- 400
burn_in_iter <- 1000
temperatures <- c(1,0.5,0.3)
bal_f <- rep("sq",length(temperatures))
#PT_IIT_sim(int p,int startsim,int endsim, int numiter,int iterswap,int burn_in, vec temp, const std::vector<std::string>& bal_function, bool bias_fix, int initial_state)

rev <- PT_IIT_sim(p,1,5,total_iter, iter_swap, burn_in_iter, temperatures, bal_f, TRUE, 500)
rev2 <- PT_a_IIT_sim(p,1,5,trunc(total_iter/iter_swap), iter_swap, burn_in_iter, temperatures, bal_f, 500, 0.0005, "iterations")
rev3 <- PT_a_IIT_sim_RF(p,1,5,total_iter, iter_swap, burn_in_iter, temperatures, bal_f, TRUE, 500, 0.0005, "iterations")


##### Testing L1_distance function #####
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
# setwd('..')
Rcpp::sourceCpp("functions/cpp_functions_highdim.cpp")

op1 <- c(0,1,1,1,0,0,0,0,0,1,1,0,0,1,0,1,1,0,0,0,1,1,0,0,0,0,1,0,1,0,1,1,1,0,0,1,0,0,1,1,0,1,1,1,1,1,1,1,0,1,0,1,0,1,0,1,0,1,0,0,1,1,0,0,0,0,0,0,1,1,0,0,1,0,1,1,1,1,0,0,0,1,1,1,0,1,1,1,0,0,0,1,0,1,1,1,1,1,0,1,1,0,1,1,0,1,0,0,1,1,0,1,1,1,0,0,0,0,1,1,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,1,0,1,1,0,0,1,1,1,1,1,0,1,1,0,1,0,1,0,1,0,1,1,0,1,0,1,1,0,0,1,0,1,0,1,1,0,1,1,1,0,0,1,1,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0,1,0,0,0,1,1,1,0,0,0,1,0,0,1,0,1,1,1,1,0,0,0,1,0,0,0,1,1,1,1,0,1,0,1,0,1,1,1,1,0,0,1,1,0,0,0,0,0,0,0,1,1,1,1,1,0,1,1,1,1,1,1,1,0,0,0,1,0,0,1,0,0,0,1,0,0,0,1,1,1,1,1,0,1,0,1,1,1,1,1,0,1,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,1,1,1,0,0,1,1,0,1,0,0,0,1,1,0,0,1,1,0,1,1,1,0,0,0,1,1,1,1,0,0,1,1,1,1,0,1,0,0,1,1,0,0,0,1,0,1,1,0,0,1,0,1,1,0,1,0,0,1,0,1,0,1,1,0,1,0,0,1,0,0,1,0,1,0,1,1,1,0,1,0,0,1,0,1,0,1,1,0,1,0,1,0,1,1,0,0,1,1,0,1,0,0,0,1,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,1,0,1,0,1,1,1,0,1,0,0,0,1,0,1,0,0,1,1,0,1,0,1,1,0,1,1,0,0,1,1,1,1,1,0,0,1,0,0,0,1,1,0,0,1,1,0,0,0,1,1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0,1,1,1,0,0,0,1,0,0,1,1,1,1,0,1,1,1,1,1,1,0,0,0,1,1,1,1,1,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1,0,1,0,0,1,1,0,1,0,1,1,1,1,0,0,0,0,0,1,0,0,0,1,0,0,1,1,0,0,0,0,1,0,1,1,1,0,1,0,1,1,0,1,0,1,1,0,0,1,1,1,0,0,0,1,1,0,0,1,1,1,1,0,0,0,1,0,1,1,0,0,1,1,0,1,1,0,1,1,0,1,1,1,0,0,0,0,0,0,0,1,1,0,1,1,1,0,1,0,0,0,0,1,0,0,1,0,0,1,1,1,0,0,0,0,1,1,1,1,0,1,1,0,1,1,1,0,1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,0,1,1,1,0,1,0,1,1,1,0,1,1,1,1,1,0,1,1,0,0,1,0,0,1,0,1,0,1,0,1,0,1,0,1,0,0,0,0,1,1,1,1,0,0,1,1,0,0,0,0,0,1,0,0,1,1,1,0,1,0,0,0,1,0,0,0,0,0,1,1,0,0,1,0,0,0,1,0,0,1,0,0,1,1,0,0,1,1,0,1,1,1,1,0,1,0,0,1,1)
op2 <- c(1,0,0,0,1,1,1,1,1,0,0,1,1,0,1,0,0,1,1,1,0,0,1,1,1,1,0,1,0,1,0,0,0,1,1,0,1,1,0,0,1,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,1,0,0,1,1,1,1,1,1,0,0,1,1,0,1,1,0,0,0,1,1,1,0,0,0,1,0,0,0,1,1,1,0,1,0,0,0,0,0,1,0,0,1,0,0,1,0,1,1,0,0,1,0,0,0,1,1,1,1,0,0,1,1,1,1,1,0,0,0,1,0,1,1,1,1,1,0,1,0,0,1,1,0,0,0,0,0,1,0,0,1,0,1,0,1,0,1,0,0,1,0,1,0,0,1,1,0,1,0,1,0,0,1,0,0,0,1,1,0,0,0,1,1,1,1,0,1,0,1,1,1,1,1,0,1,0,1,1,1,0,0,0,1,1,1,0,1,1,0,1,0,0,0,0,1,1,1,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,0,1,1,1,0,1,1,1,1,1,1,1,0,0,0,0,0,1,0,0,0,1,0,0,0,1,1,1,0,1,1,0,1,1,1,0,1,1,1,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0,1,1,1,0,1,1,1,0,1,0,1,1,1,1,1,1,0,0,0,1,1,0,0,1,0,1,1,1,0,0,1,1,0,0,1,0,0,0,1,1,1,0,0,0,0,1,1,0,0,0,0,1,0,1,1,0,0,1,1,1,0,1,0,0,1,1,0,1,0,0,1,1,1,1,0,1,0,1,0,0,1,0,1,1,0,1,1,0,1,0,1,0,0,0,1,0,1,1,0,1,0,1,0,0,1,0,1,0,1,0,0,1,1,0,0,1,0,1,1,1,0,1,1,1,0,0,0,1,1,1,1,1,1,0,0,0,0,1,0,1,0,0,0,1,0,1,1,1,0,1,0,1,1,0,0,1,0,1,0,0,1,0,0,1,1,0,0,0,0,0,1,1,0,1,1,1,0,0,1,1,0,1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,1,1,1,0,0,1,1,1,0,0,1,1,1,1,0,1,1,0,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,1,1,0,1,1,1,1,1,0,1,0,1,0,1,1,0,0,1,0,1,0,0,0,0,1,1,0,1,1,0,1,1,1,0,1,1,0,0,1,1,1,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,1,1,0,0,0,1,1,1,0,0,1,1,0,0,0,0,1,1,1,0,1,0,0,1,1,0,0,1,0,0,1,0,0,1,0,0,0,1,1,1,1,1,1,1,0,0,1,0,0,0,1,0,1,1,1,1,0,1,1,0,1,1,0,0,0,1,1,1,1,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,0,0,1,0,1,1,1,1,1,1,0,1,1,1,1,1,0,0,0,1,1,1,0,1,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,1,0,0,1,1,0,1,1,0,1,1,1,0,1,0,1,0,1,0,1,1,1,1,0,0,0,0,1,1,0,0,1,1,1,1,1,0,1,1,0,0,0,1,0,1,1,1,0,1,1,1,1,1,0,0,1,1,0,1,1,1,0,1,1,0,1,0,0,0,1,1,0,0,1,0,0,0,0,1,0,1,1,0,0)

L1_distance(op1,op2)
L1_distance(c(0,1,1,1,1,0),c(0,0,1,0,0,1))

p <- 800
total_iter <- 2000
iter_swap <- 2000
burn_in_iter <- 10
temperatures <- c(1,0.5,0.1,0.05)
bal_f <- rep("sq",length(temperatures))

rev <- PT_IIT_sim(p,1,2,total_iter, iter_swap, burn_in_iter, temperatures, bal_f, TRUE, filename="gset/G1.txt", num_states_visited=5, starting_coord=-1)
rev2 <- PT_a_IIT_sim_RF(p,1,2,total_iter, iter_swap, burn_in_iter, temperatures, bal_f, TRUE, filename="gset/G1.txt", num_states_visited=5, starting_coord=-1, decreasing_constant=.5, reduc_model="never")
rev3 <- PT_a_IIT_sim(p,1,1,trunc(total_iter/iter_swap), iter_swap, burn_in_iter, temperatures, bal_f,filename="gset/G1.txt", num_states_visited=5, starting_coord=-1, decreasing_constant=.5, reduc_model="never")


p <- 800
total_iter <- 200
iter_swap <- total_iter
burn_in_iter <- 1
temperatures <- c(.1)
# temperatures <- c(0)
bal_f <- rep("sq",length(temperatures))
rev3 <- PT_a_IIT_sim(p,1,1,trunc(total_iter/iter_swap), iter_swap, burn_in_iter, temperatures, bal_f,filename="gset/G1.txt", num_states_visited=5, starting_coord=-1, decreasing_constant=0, reduc_model="never")
m1 <- rev3$distance_mode1[,1,1]
m2 <- rev3$distance_mode2[,1,1]
m3 <- rev3$distance_origin[,1,1]

tibble(y1=m1,y2=m2,y0=m3,x=1:length(m1)) |> 
  pivot_longer(y1:y0,names_to = "mode",values_to = "distance") |> 
  ggplot(aes(x=x,y=distance, col=mode))+
  geom_line()+labs(title='A-IIT')

rev <- PT_IIT_sim(p,1,1,total_iter, iter_swap, burn_in_iter, temperatures, bal_f, TRUE, filename="gset/G1.txt", num_states_visited=5, starting_coord=-1)
m1 <- rev$distance_mode1[,1,1]
m2 <- rev$distance_mode2[,1,1]
m3 <- rev$distance_origin[,1,1]

tibble(y1=m1,y2=m2,y0=m3,x=1:length(m1)) |> 
  pivot_longer(y1:y0,names_to = "mode",values_to = "distance") |> 
  ggplot(aes(x=x,y=distance, col=mode))+
  geom_line()+labs(title='PT-IIT')

rev2 <- PT_a_IIT_sim_RF(p,1,1,total_iter, iter_swap, burn_in_iter, temperatures, bal_f, TRUE, filename="gset/G1.txt", num_states_visited=5, starting_coord=-1, decreasing_constant=.5, reduc_model="never")
m1 <- rev2$distance_mode1[,1,1]
m2 <- rev2$distance_mode2[,1,1]
m3 <- rev2$distance_origin[,1,1]

tibble(y1=m1,y2=m2,y0=m3,x=1:length(m1)) |> 
  pivot_longer(y1:y0,names_to = "mode",values_to = "distance") |> 
  ggplot(aes(x=x,y=distance, col=mode))+
  geom_line()+labs(title='PT-IIT RF')

#### Testing L1_distance function in bimodal highdim example
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("functions/cpp_functions_highdim_2.cpp")
p <- 200
numiter <- 5000
burn_in <- 1000
iterswap <- 5000
total_swaps <- numiter/iterswap
Q_matrix <- matrix(0,nrow=p,ncol=2)
for(i in 1:p){
  if(i%%2==1){Q_matrix[i,2]=1}
  if(i%%2==0){Q_matrix[i,1]=1}
}
mode1 <- Q_matrix[,1]
mode2 <- Q_matrix[,2]
L1_distance(mode1,mode2)
L1_distance(mode1,rep(0,p))

temperatures <- c(0.1)
bal_f <- rep("sq",length(temperatures))

set.seed(456)
check_res <- PT_IIT_sim(p,1,1,numiter, iterswap,burn_in, temperatures,bal_f, TRUE,"anything",3,-1)
output <- check_res
tibble(y1=output$distance_mode1,y2=output$distance_mode2,y0=output$distance_origin,x=1:length(output$distance_mode1)) |> 
  pivot_longer(y1:y0,names_to = "mode",values_to = "distance") |> 
  ggplot(aes(x=x,y=distance, col=mode))+
  geom_line()+labs(title='PT-IIT RF')

set.seed(456)
check_ad <- PT_a_IIT_sim(p,1,1,total_swaps,iterswap,burn_in, temperatures, bal_f,"anything",3,-1, decreasing_constant = 0, reduc_model = "never")
output <- check_ad
tibble(y1=output$distance_mode1,y2=output$distance_mode2,y0=output$distance_origin,x=1:length(output$distance_mode1)) |> 
  pivot_longer(y1:y0,names_to = "mode",values_to = "distance") |> 
  ggplot(aes(x=x,y=distance, col=mode))+
  geom_line()+labs(title='PT-IIT RF')

set.seed(456)
check_rf <- PT_a_IIT_sim_RF(p,1,1,numiter,iterswap,burn_in, temperatures, bal_f,TRUE, "anything",3,-1, decreasing_constant = 0, reduc_model = "never")
output <- check_rf
tibble(y1=output$distance_mode1,y2=output$distance_mode2,y0=output$distance_origin,x=1:length(output$distance_mode1)) |> 
  pivot_longer(y1:y0,names_to = "mode",values_to = "distance") |> 
  ggplot(aes(x=x,y=distance, col=mode))+
  geom_line()+labs(title='PT-IIT RF')


### Checking the possible bias of A-IIT
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("functions/cpp_functions_lowdim_2.cpp")
p <- 16
Q_matrix <- matrix(0,nrow=p,ncol=2)
for(i in 1:p){
  if(i%%2==0){Q_matrix[i,2]=1}
  if(i%%2==1){Q_matrix[i,1]=1}
}
mode1 <- Q_matrix[,1]
mode2 <- Q_matrix[,2]


a <- a_IIT_update(mode1,"sq",1,0,FALSE, 0, 0,0)
b <- a_IIT_update(mode2,"sq",1,0,FALSE, 0, 0,0)

loglik(mode1)
l1 <- c()
l2 <- c()
for(i in 1:p){
  tempx <- mode1;
  tempx[i] <- 1-tempx[i];
  l1[i] <- loglik(tempx)
  tempx <- mode2;
  tempx[i] <- 1-tempx[i];
  l2[i]<- loglik(tempx)
}

#For this parameters the normalizing constant is
((1+exp(-6))^p + (1+exp(-6))^p)

#Con estos calculos, usar como bounding constant gamma=0
#nos da lo mismo que cambiar la balancing function
# a ser MIN{1,r}, o sea es traditional MH
#El problema era que cuando se me acababan los L samples
#Yo si cambiaba de estado, y no debia


##### Testing function to find temperatures #####
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("functions/find_temp_func.cpp")

p <- 16
interswap <- 50
theta <- 6
temp_ini <- 1
bal_function <- "sq"

set.seed(123)
temperature_PT_IIT(p,interswap, temp_ini,bal_function, theta)



p <- 800
interswap <- 100
temp_ini <- 1
bal_function <- "sq"
theta <- 3

set.seed(125)
temperature_PT_IIT(p,interswap, temp_ini,bal_function, theta)
# FINAL RESULTS:
#   Swap: 374 avg. swap prob: 0.234261 new temperature: 0.90526


set.seed(888)
temperature_PT_IIT(p,interswap, temp_ini,bal_function, theta)


rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("functions/find_temp_func.cpp")

p <- 16
interswap <- 100
theta <- 6
temp_ini <- 1
bal_function <- "sq"

set.seed(13)
temperature_PT_a_IIT(p,interswap, temp_ini,bal_function, theta)
# FINAL RESULTS:
#   Swap: 97553 avg. swap prob: 0.234799 new temperature: 0.380702


set.seed(54)
temperature_PT_a_IIT(p,interswap, temp_ini,bal_function, theta)
# FINAL RESULTS:
#   Swap: 22866 avg. swap prob: 0.234593 new temperature: 0.082487

#log(1/x - 1)
#1/(1+exp(x))


##### Testing writing files

for(i in 1:5){
  if(i==1){
    write(paste0("Iteration \n",i," result: ",runif(1)), file = "results/results.txt", append = FALSE)
  }else{
    write(paste0("Iteration ",i," result: ",runif(1)), file = "results/results.txt", append = TRUE)
  }
}

##### Checking number of threads #####

library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("functions_other/testing_cpp_functions.cpp")

check_threads()

check_rng()


#### Comparing Parallelization of IIT_update ###
Rcpp::sourceCpp("functions_other/testing_cpp_functions.cpp")

library(rbenchmark)
p <- 2000
Q_matrix <- matrix(0,nrow=p,ncol=2)
for(i in 1:p){
  if(i%%2==0){Q_matrix[i,2]=1}
  if(i%%2==1){Q_matrix[i,1]=1}
}
chosen_bf <- "sq"
temperature <- 1
theta <- 3
seed <- 625
X <- rbinom(p,1,0.5)


IIT_update_w(X,Q_matrix,chosen_bf,temperature,theta,seed)


# Then test with parallel
library(OpenMPController)
omp_set_num_threads(12)  # Use 4 threads
asd <- IIT_update_w_parallel(X,Q_matrix,chosen_bf,temperature,theta,seed)

res <- benchmark(IIT_update_w_parallel(X,Q_matrix,chosen_bf,temperature,theta,seed),IIT_update_w(X,Q_matrix,chosen_bf,temperature,theta,seed),replications = 500,order="relative")


######### Test for each parallel ############
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
# Sys.setenv("PKG_CXXFLAGS" = "-std=c++17 -DHAS_PARALLEL")
# Sys.setenv("CXX_STD" = "CXX17")

Rcpp::sourceCpp("functions_other/test_parallel.cpp", verbose = TRUE)

# Rcpp::sourceCpp("functions/find_temp_func.cpp", verbose = TRUE)


# run_sim(int p, const int iter,const vec temp,const std::string bal_function,
#         double theta=3,int seed=5)
p <- 16
iter <- 50
temp <- c(1,0.5)
theta <- 3
seeed <- 123
set.seed(5)
X <- run_sim(p,iter,temp,"sq",theta,seeed)
set.seed(6)
X2 <- run_sim(p,iter,temp,"sq",theta,seeed)


######### Test RcppParallel ############
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
Rcpp::sourceCpp("functions/find_temp_parallel.cpp", verbose = TRUE)

p <- 500
num_iter <- 500
temperatures <- c(1,0.8,0.2,0.3,0.5)
bal_func <- "sq"
theta <- 3

set.seed(123)
c <- PT_IIT_parallel_sim(p,num_iter,temperatures,bal_func,theta)
c2 <- PT_IIT_parallel_sim(p,num_iter,temperatures,bal_func,theta)
identical(c,c2)
set.seed(123)
d <- PT_IIT_parallel_sim_wfor(p,num_iter,temperatures,bal_func,theta)

library(rbenchmark)

results <- benchmark(PT_IIT_parallel_sim(p,num_iter,temperatures,bal_func,theta),
                     PT_IIT_parallel_sim_wfor(p,num_iter,temperatures,bal_func,theta),
                     replications=200)


#### Testing little by little
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
Rcpp::sourceCpp("functions/find_temp_parallel.cpp", verbose = TRUE)
p <- 20
temperature <- 1
bal_func <- 2;
theta <- 3;
set.seed(30)
X <- test_1(p,temperature,bal_func,theta)
Y <- test_1(p,temperature,bal_func,theta)
set.seed(30)
Y <- test_1(p,temperature,bal_func,theta)
identical(X,Y)
#sum(exp(X[1:(length(X)-1)]))

############## Testing find_temps function
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
Rcpp::sourceCpp("functions/find_temp_parallel.cpp")

#temperature_PT_IIT(int p,int interswap, double temp_ini, int bal_func, const double& theta)
p <- 100
interswap <- 10
temp_ini <- 1
bal_func <- 2
theta <- 3
set.seed(45)
result1 <- temperature_PT_IIT(p,interswap,temp_ini,bal_func,theta)


################### Checking el código de Marco
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
library(RcppProgress)
Rcpp::sourceCpp("C:/Users/ralex/Downloads/toy-example-target-mixture-mv.cpp")


# pt_cmh_parallel(
#   int nsim,
#   NumericVector init, 
#   NumericVector temp_vector,
#   NumericVector means,
#   NumericVector weights,
#   arma::Col<double> step_size,
#   NumericVector sd,
#   bool step_depends_temp,
#   int within_temp,
#   bool display_progress=false)



test <- pt_cmh_parallel(100000,
                c(1,1),
                c(1,0.8),
                c(1,2),
                c(0.5,0.5),
                c(0.2,0.2),
                c(0.5,0.5),
                TRUE,
                500100,
                T)

########### Testing Gibbs Sampler
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
Rcpp::sourceCpp("functions/gibbs_sampler_temp.cpp")





# run_gibbs(int p, int num_iter, int burn_in, double temp_ini,double theta, int base_seed, bool adapting_factors)
p <- 10
num_iter <- 1
burn_in <- 100
temp_ini <- 1
bal_func <- 2
theta <- 3
base_seed <- 123
adapting_factors <- TRUE
set.seed(base_seed)
results <- find_temp_gibbs(p,num_iter,burn_in,temp_ini,bal_func,theta,base_seed,adapting_factors)

find_temp_gibbs(p,num_iter,burn_in,temp_ini,bal_func,theta,base_seed,F)


test_flip_coord(1,F,3)


########### Testing NEW Gibbs Sampler
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
Rcpp::sourceCpp("functions/find_temp_parallel.cpp", verbose=T)




p <- 1000
interswap <- 1
burn_in <- 5000
temp_ini <- 1
bal_func <- 2
theta <- 3
base_seed <- 123
adapting_factors <- TRUE
set.seed(base_seed)
# (int p,int interswap, double temp_ini, int bal_func, const double& theta, int base_seed)
results <- find_temp_gibbs(p,interswap,burn_in,temp_ini,bal_func,theta,base_seed,adapting_factors)
set.seed(base_seed)
results2 <- find_temp_gibbs(p,interswap,burn_in,temp_ini,bal_func,theta,base_seed,adapting_factors)

rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
Rcpp::sourceCpp("functions/find_temp_parallel.cpp")
# find_temp_gibbs_PT_IIT(int p, int burn_in,double temp_ini, int bal_func, const double& theta, int gibbs_steps)
p <- 2000
burn_in <- 10
temp_ini <- 1
bal_func <- 2
theta <- 3
base_seed <- 123
gibbs_steps <- 100
set.seed(45)# iF rho=0 converges fast,p=16
set.seed(554)
results <- find_temp_gibbs_PT_IIT(p,burn_in,temp_ini,bal_func,theta,gibbs_steps,1)

interswap <- 1
p <- 1000
set.seed(45)# iF rho=0 converges fast, p=16
results_a_iit <- find_temp_gibbs_A_IIT(p,interswap,burn_in,temp_ini,bal_func,theta,base_seed)


### Testing parallelized algorithm
rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
Rcpp::sourceCpp("functions/highdim_2_parallel.cpp")

# PT_IIT_sim(int p,int startsim,int endsim, int numiter, int iterswap,int burn_in, vec temp, int bal_func, bool bias_fix,const std::string& filename,int num_states_visited,const std::vector<int>& starting_coord, double theta)
p <- 1000
startsim <- 1
endsim <- 5
numiter <- 40000
interswap <- 100
burn_in <- 10000
temp <- c(1,0.917675,0.846553,0.78432,0.725)
# temp <- c(1,0.918)
bal_func <- 2
bias_fix <- T
file_name <- ""
states_visit <- 2
starting_coord <- c(0.0)
theta <- 3
set.seed(90)
results <- PT_IIT_sim(p,startsim,endsim,numiter,interswap,burn_in,temp,
           bal_func,bias_fix,file_name,states_visit,starting_coord,theta)

Q_matrix <- matrix(0,nrow=p,ncol=2)
for(i in 1:p){
  if(i%%2==1){Q_matrix[i,2]=1}
  if(i%%2==0){Q_matrix[i,1]=1}
}

loglik(Q_matrix[,1],Q_matrix,3)
loglik(c(1,Q_matrix[2:p,1]),Q_matrix,3)
loglik(c(rep(0,p-1),1),Q_matrix,theta)


p <- 1000
startsim <- 1
endsim <- 3
numiter <- 200000
interswap <- 50
burn_in <- 2000
temp <- c(1.639838,1.472027,1.327346,1.202306,1.194742,1.092,1)
# temp <- c(1,0.918)
bal_func <- 2
bias_fix <- T
file_name <- ""
states_visit <- 2
starting_coord <- c(0.0)
theta <- 3
set.seed(runif(1))
results <- PT_IIT_sim(p,startsim,endsim,numiter,interswap,burn_in,temp,
                      bal_func,bias_fix,file_name,states_visit,starting_coord,theta)

results$distance_mode1
results$distance_mode2
results$time_mode1
results$time_mode2
results$time_taken
results$swap_rate



set.seed(123)
results_123 <- PT_IIT_sim(p,startsim,endsim,numiter,interswap,burn_in,temp,
                      bal_func,bias_fix,file_name,states_visit,starting_coord,theta)
results_bk <- results_123;
set.seed(123)
results <- PT_IIT_sim(p,startsim,endsim,numiter,interswap,burn_in,temp,
                          bal_func,bias_fix,file_name,states_visit,starting_coord,theta)
#compare results with results_123
results$swap_rate;results_123$swap_rate;
results$distance_mode1; results_123$distance_mode1;
results$distance_mode2; results_123$distance_mode2;


set.seed(123123)
results_123123 <- PT_IIT_sim(p,startsim,endsim,numiter,interswap,burn_in,temp,
                      bal_func,bias_fix,file_name,states_visit,starting_coord,theta)

Rcpp::sourceCpp("functions_other/testing_cpp_functions.cpp")
set.seed(43)
(vv <- sample(1:10,replace=F))
check_min_find(vv,4)
