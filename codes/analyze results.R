### Process simulation results
library(stringr)
library(tidyverse)
# data_sum <- data.frame(file=dir("results"))
# data_sum$alg <- str_extract(results, ".*(?=_sim)")


data_sum <- tibble(file = dir("results")) |> 
  mutate(alg = str_extract(file, ".*(?=_sim)"),
         sim = str_extract(file,"(?<=sim_)[^_]+"),
         iter=str_extract(file,"(?<=iter_)[^(_\\.)]+"),
         interswap=str_extract(file,"(?<=interswap_)[^_]+"),
         tot_swap=str_extract(file,"(?<=totalswap_)[^\\.]+"),
         iterswap=str_extract(file,"(?<=iterswap_)[^\\.]+")) |> 
  mutate(across(-(1:2), as.numeric))


for(i in 1:nrow(data_sum)){
  data <- readRDS(file.path("results",data_sum[i,1]))
  print(data_sum[i,"alg"])
  print(names(data))
}

