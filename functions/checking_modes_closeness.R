library(tidyverse)
num_boxes <- 2
boxes <- 1:num_boxes
B <- 50000
max_balls <- 10
vec_count <- c()
skip <- 3

for(balls in (num_boxes+skip):max_balls){
  print(paste0("Balls = ",balls))
  counting <- 0
  for(i in 1:B){
    checa <- sample(boxes,size=balls,replace=TRUE)
    if(length(unique(checa))==num_boxes){
      counting <- counting+1;
    }
  }
  vec_count[balls] <- counting/B
}
vec_count <- vec_count[(num_boxes+skip):max_balls]

data <- tibble(balls=(num_boxes+skip):max_balls,prob=vec_count)

data |> ggplot(aes(x=balls,y=prob))+
  geom_point()
tail(data,17)
# plot(x=(num_boxes+2):20,y=vec_count)
# counting/B

### For 7 modes, 21 balls give +0.74 of probability
### For 7 modes, 30 balls give +0.92 of probability
### For 7 modes, 35 balls give +0.95 of probability
### For 7 modes, 45 balls give +0.99 of probability

### For 5 modes, 19 balls give +0.90 of probability
### For 5 modes, 21 balls give +0.95 of probability
### For 5 modes, 25 balls give +0.98 of probability
### For 5 modes, 29 balls give +0.99 of probability

### For 2 modes, 5 balls give +0.93 of probability
### For 2 modes, 6 balls give +0.95 of probability
### For 2 modes, 10 balls give +0.99 of probability
