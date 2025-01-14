##### General functions #####


#Transform nummbers in base 10 to a vector representation of the number in base 2
#Inputs:
### Number to transform
### Maximum position needed in base 2 (dimension of the vector)
#Output: Vector of 0-1 representing the number in base 2
NumberToVector <- function(number,N)
{
  if(number < 0 || number >= 2^N)
  {
    print("Error, number bigger than dimension")
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

#Computes Total variation distance
#Inputs:
### Vectors with density to compare
# It standardize the vectors if needed
#Output:
### A value between 0 and 1 indicating the TVD
TVD <- function(pi,pi.est){
  if(sum(pi)!=1){pi <- pi/sum(pi)}
  if(sum(pi.est)!=1){pi.est <- pi.est/sum(pi.est)}
  return(0.5*sum(abs(pi-pi.est)))
}





##### Functions to transform output from PT simulation #####

#Compute roundtrip rate
#Inputs:
### Index process (from Rcpp function)
### Number oof iterations, total simulations and iteration between swaps.
#Output: matrix reporting the number of Round Trips for each simulation
PT_RT <- function(ip,total_iter,iterswap, total_sim){
  ip <- as.data.frame(ip)
  total_swaps <- floor(total_iter/iterswap)
  if(nrow(ip)!=(total_sim*total_swaps+1)){print("Error with size of input")}
  
  max_temp <- max(ip[1,])#Get max temp
  min_temp <- min(ip[1,])#Get min temp
  
  initial_ip <- ip[1,]#Get the initial state of the index process

  ip <- ip |> slice(-1)
  ip$sim <- rep(1:total_sim,each=total_swaps) #Add a column labeling each row with the simulation it corresponds
  ip <- ip |> select(sim,everything())# Re-order columns so simulation number is first column
  
  #The round trip starts with the hottest temperature
  
  rt_count <- matrix(0,nrow = total_sim,ncol=ncol(ip)-1) # to track round trip rate
  
  for (c in 2:ncol(ip)){#Loop for each replica, ignoring first column that contains the number of simulation
    first_last <- (initial_ip[c-1]==max_temp)*max_temp + (initial_ip[c-1]==min_temp)*min_temp
    for(i in 1:nrow(ip)){#Loop for each iteration
      if(i==1){# Restart the index process for each simulation
        last <- first_last
        start_rt <- FALSE
      }else if(ip[i,1]!=ip[i-1,1]){#If we changed the number of simulation
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
  # Code to test  
  # ip |> filter(V4==1|V4==4) |> select(sim,V4) |> filter(sim==1)
  # ip |> filter(V1==1|V1==4) |> select(sim,V1) |> filter(sim==1)
  # ip |> filter(V2==1|V2==4) |> select(sim,V2) |> filter(sim==1)
  # ip |> filter(V3==1|V3==4) |> select(sim,V3) |> filter(sim==1)
  # ip |> filter(V3==1|V3==4) |> select(sim,V3) |> filter(sim==9)
  # ip |> filter(V1==1|V1==4) |> select(sim,V1) |> filter(sim==9)
  # ip|> filter(sim==9) |> pull(V1)
  # ip |> filter(sim==9) |> pull(V3)
  return(rt_count);
}



#Compute likelihood based on distance to a list of modes
#Inputs:
#Vector X to compute the likelihood
#List of modes (same dimension as the vector)
#Theta parameter to multiply th exponent -> bigger theta implies more concentration of probability around modes
#dimension of the vector p
ll_comp <- function(X,modes_list, theta=1,p){
  if(length(X)!=p){print('Error in dimension');return(-1);}
  total=0;
  for(m in 1:length(modes_list)){# loop for modes
    total = total+exp(-theta*sum(abs(X-modes_list[[m]])))
  }
  return(log(total));
}