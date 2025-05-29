sq_b <- function(x,bb){
  boundf <- min(sqrt(x),bb)/bb
  boundf_1 <- min(sqrt(1/x),bb)/bb
  
  return(min(boundf,x*boundf_1))
}

TVD <- function(pi,pi.est){
  return(sum(abs(pi-pi.est))*0.5)
}
# Defining target distribution pi
space <- 1:10
pi.dist <- rep(1,length(space))
pi.dist[2] <- 2
pi.dist[7] <- 7

pi.true <- pi.dist/sum(pi.dist)

# Parameters
bound <- 0.1
pi.est <- rep(0,length(space))

tvd <- 1 #initialize TVD
iter.count <- 0
x.current <- sample(1:10,size=1) #Initialize X
while(tvd>0.2){
  iter.count <- iter.count+1
  if(x.current==1){ #Only neighbor is 2
    Z_h <- sq_b(pi.dist[2]/pi.dist[x.current],bound)
    pi.est[x.current] <- pi.est[x.current]+1+rgeom(1,Z_h)
    x.current <- 2
  }else   if(x.current==10){ #Only neighbor is 9
    Z_h <- sq_b(pi.dist[9]/pi.dist[x.current],bound)
    pi.est[x.current] <- pi.est[x.current]+1+rgeom(1,Z_h)
    x.current <- 9
  } else{ #Consider 2 neighbor states
    Z_h <- mean(sq_b(pi.dist[x.current-1]/pi.dist[x.current],bound),sq_b(pi.dist[x.current+1]/pi.dist[x.current],bound))
    pi.est[x.current] <- pi.est[x.current]+1+rgeom(1,Z_h)
    x.current <- sample(c(x.current-1,x.current+1), size=1)
  }
  
  tvd <- TVD(pi.true,pi.est/sum(pi.est))
  print(tvd)
}
iter.count
pi.est/sum(pi.est)
pi.true

######## Plot
x <- c(1/7,1/2,1,2,7)
data_u <- tibble(r=x,
                 `9`=sapply(x,sq_b,bb=9),
                 `1`=sapply(x,sq_b,bb=1),
                 `2`=sapply(x,sq_b,bb=2),
                 `3`=sapply(x,sq_b,bb=3),
                 `4`=sapply(x,sq_b,bb=4),
                 `5`=sapply(x,sq_b,bb=5))

data_u |> pivot_longer(-r,names_to = 'b.fun',values_to = 'h(r)') |>  
  ggplot(aes(x=r,y=`h(r)`,color = `b.fun`))+
  geom_line(size=0.2)+
  geom_point()+
  labs(col=TeX("$\\gamma$"),
       y=TeX("$h_{\\gamma}(r)$"),
       x=TeX("$r$"),
       title=TeX("$\\sqrt{r}$ with different bounds"))


######## Computing Z_h factors


bound <- 2.5
for(i in 1:10){
  x.current <- i
  if(x.current==1){ #Only neighbor is 2
    Z_h <- sq_b(pi.dist[2]/pi.dist[x.current],bound)
  }else   if(x.current==10){ #Only neighbor is 9
    Z_h <- sq_b(pi.dist[9]/pi.dist[x.current],bound)
  } else{ #Consider 2 neighbor states
    Z_h <- mean(sq_b(pi.dist[x.current-1]/pi.dist[x.current],bound),sq_b(pi.dist[x.current+1]/pi.dist[x.current],bound))
  }
  print(paste0(Z_h,", ",1/Z_h))
}

######## Traditional MH for benchmark

# Parameters
pi.est <- rep(0,length(space))

tvd <- 1 #initialize TVD
iter.count <- 0
x.current <- sample(1:10,size=1) #Initialize X
while(tvd>0.06){
  iter.count <- iter.count+1
  if(x.current==1){ #Only neighbor is 2
    Z_h <- min(pi.dist[2]/pi.dist[x.current],1)
    pi.est[x.current] <- pi.est[x.current]+1+rgeom(1,Z_h)
    x.current <- 2
  }else   if(x.current==10){ #Only neighbor is 9
    Z_h <- min(pi.dist[9]/pi.dist[x.current],1)
    pi.est[x.current] <- pi.est[x.current]+1+rgeom(1,Z_h)
    x.current <- 9
  } else{ #Consider 2 neighbor states
    Z_h <- mean(min(pi.dist[x.current-1]/pi.dist[x.current],1),min(pi.dist[x.current+1]/pi.dist[x.current],1))
    pi.est[x.current] <- pi.est[x.current]+1+rgeom(1,Z_h)
    x.current <- sample(c(x.current-1,x.current+1),size=1)
  }
  
  tvd <- TVD(pi.true,pi.est/sum(pi.est))
  # print(tvd)
}
