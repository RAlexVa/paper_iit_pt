library(tidyverse)
library(latex2exp)
library(here)
###### Bounded balancing functions ######
bf <- function(x,f,bb){
  boundf <- min(f(x),bb)/bb
  boundf_1 <- min(f(1/x),bb)/bb
  
  return(min(boundf,x*boundf_1))
}

bf_nobb <- function(x,f,bb){
  boundf <- min(f(x),bb)
  boundf_1 <- min(f(1/x),bb)
  
  return(min(boundf,x*boundf_1))
}

# for(i in seq(0,15,by=0.1)){
#   print(paste(round(bf(i,sqrt,4),4),round(bf(1/i,sqrt,4),4),round(i*bf(1/i,sqrt,4),10)==round(bf(i,sqrt,4),10)))
# }

x <- seq(0.1,20,by=0.1)
#### Random balancing functions
data_u <- tibble(r=x,
                 f1=sapply(x,bf,f=sqrt,bb=4),
                 f2=sapply(x,bf,f=sqrt,bb=2),
                 f3=sapply(x,bf,f=function(x){x^(2)},bb=2),
                 f4=sapply(x,bf,f=function(x){x^(1/4)},bb=2),
                 f5=sapply(x/16,min,1),
                 f6=sapply(x/5,min,1))
data_u |> pivot_longer(-r,names_to = "b.fun",values_to = 'h(r)') |>  
  ggplot(aes(x=r,y=`h(r)`,color = `b.fun`))+
  geom_line(size=1)




x <- seq(0.1,5,by=0.1)
#### Random balancing functions
data_u <- tibble(r=x,
                 f1=sapply(x,bf,f=sqrt,bb=1.5),
                 f2=sapply(x,bf,f=sqrt,bb=2))
data_u |> pivot_longer(-r,names_to = "b.fun",values_to = 'h(r)') |>  
  ggplot(aes(x=r,y=`h(r)`,color = `b.fun`))+
  geom_line(size=1)


# data_u <- tibble(r=x,
#                  f1=sapply(x,bf,f=function(x){1+x},bb=.1),
#                  f2=sapply(x,bf,f=function(x){1+x},bb=1),
#                  f3=sapply(x,bf,f=function(x){1+x},bb=2),
#                  f4=sapply(x,bf,f=function(x){1+x},bb=3),
#                  f5=sapply(x,bf,f=function(x){1+x},bb=4),
#                  f6=sapply(x,bf,f=function(x){1+x},bb=5))
# 
# data_u |> pivot_longer(-r,names_to = "b.fun",values_to = 'h(r)') |>  
#   ggplot(aes(x=r,y=`h(r)`,color = `b.fun`))+
#   geom_line(size=1)


data_u <- tibble(r=x,
                 `1`=sapply(x,bf,f=sqrt,bb=1),
                 `2`=sapply(x,bf,f=sqrt,bb=2),
                 `3`=sapply(x,bf,f=sqrt,bb=3),
                 `4`=sapply(x,bf,f=sqrt,bb=4),
                 `5`=sapply(x,bf,f=sqrt,bb=5),
                 `6`=sapply(x,bf,f=sqrt,bb=6))

data_u |> pivot_longer(-r,names_to = 'b.fun',values_to = 'h(r)') |>  
  ggplot(aes(x=r,y=`h(r)`,color = `b.fun`))+
  geom_line(size=1)+
  labs(col=TeX("$\\gamma$"),
       y=TeX("$h_{\\gamma}(r)$"),
       x=TeX("$r$"),
       title=TeX("$\\sqrt{r}$ with different bounds"))


data_u <- tibble(r=x,
                 f1=sapply(x,bf,f=sqrt,bb=4),
                 f2=sapply(x/4^2,min,1),
                 f3=sapply(x/16,bf_nobb,f=function(x){x^2},bb=4))

data_u |> pivot_longer(-r,names_to = 'b.fun',values_to = 'h(r)') |>  
  ggplot(aes(x=r,y=`h(r)`,color = `b.fun`))+
  geom_line(size=1)+
  labs(col=TeX("b.fun"),
       y=TeX("$h(r)$"),
       x=TeX("$r$"),
       title=TeX("Possible behaviors of bounded balancing functions"))


# Using a power function doesnt seems to give good results
#But perhaps we can just use min function with a different bounding constant
#Here I think we can choose how important we want to make each neighbors


###### unbounded balancing functions ######
x <- seq(0,4,by=0.1)
data_u <- tibble(r=x,
                 min=sapply(x,min,1),
                 sq=sqrt(x),
                 max=sapply(x,max,1))

(plot1 <- data_u |> pivot_longer(-r,names_to = "h",values_to = 'h(r)') |>  
    ggplot(aes(x=r,y=`h(r)`,color = `h`))+
    geom_line(size=1)+
    geom_segment(aes(x=2,y=0,xend=2,yend=2),color = "blue", linetype = "dashed", linewidth = 1)+
    geom_segment(aes(x=3,y=0,xend=3,yend=3),color = "red", linetype = "dashed", linewidth = 1)+
    annotate("text", x=2.1, y=0, label= TeX("$y_1$"),size=6, color='blue')+
    annotate("text", x=3.1, y=0, label= TeX("$y_2$"),size=6, color='red')+
    theme_minimal()+
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=18),
          legend.title = element_text(size = 25),
          legend.text=element_text(size = 20)))

jpeg(file.path(here(),"fig","compare_balancing.jpg"), width = 850, height = 300)
plot1
dev.off()

###### Bounded version of sqrt balancing function ######
x <- seq(0.1,6,by=0.1)
#### Random balancing functions
data_u <- tibble(r=x,
                 `1`=sapply(x,bf,f=sqrt,bb=1),
                 `1.3`=sapply(x,bf,f=sqrt,bb=1.3),
                 `1.5`=sapply(x,bf,f=sqrt,bb=1.5),
                 `1.7`=sapply(x,bf,f=sqrt,bb=1.7),
                 `2`=sapply(x,bf,f=sqrt,bb=2),
                 `2.5`=sapply(x,bf,f=sqrt,bb=2.5),
                 `3`=sapply(x,bf,f=sqrt,bb=3))
(plot_sqrt <- data_u |> pivot_longer(-r,names_to = "b.fun",values_to = 'h(r)') |>  
  ggplot(aes(x=r,y=`h(r)`,color = `b.fun`))+
  geom_line(size=1)+
  labs(y = TeX("$h_{\\gamma}(r)$"),color=TeX("$\\gamma$"))+
  theme_minimal()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18),
        legend.title = element_text(size = 25),
        legend.text=element_text(size = 20)))

jpeg(file.path(here(),"fig","bounded_sqrt.jpg"), width = 850, height = 300)
plot_sqrt
dev.off()
