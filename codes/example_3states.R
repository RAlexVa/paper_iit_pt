##### Simple example of RF algorithm #####
set.seed(2412)
states <- c(1,2,3) #state space
freq <- c(0,0,0) #To store how many times each state is visited
N <- 3*10^5 #Total number of samples
x <- sample(1:3,size=1) #Initialize x
for(i in 1:N){
  state <- sample(states[-x],size=1) #sample uniformly other states
  freq[state] <- freq[state]+1 #Update frequency
  x <- state #update state
}

# Estimated density
freq/N

#Weighted density
freq*c(1,2,1)/sum(freq*c(1,2,1))
