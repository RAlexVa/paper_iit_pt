library(tidyverse)
# setwd("C:/Users/ralex/Documents/src/paper_iit_pt")
source("codes/lowd_seeded.R")

list_ids <- 528:587
# id_chosen <- 400;
swap_rate <- tibble(id=numeric(),r1=numeric(),r2=numeric(),r3=numeric())
iterations <- tibble(id=numeric(),r1=numeric(),r2=numeric(),r3=numeric(),r4=numeric())
for(id_chosen in list_ids){
  print(paste0("Running id:",id_chosen))
  run_lowd(id_chosen)
  
  lll <- readRDS(file.path(getwd(),"results",paste0("sim_lowdim_id_",id_chosen,"_1.Rds")))
  
  swap_rate_temp <- tibble(
    id=id_chosen,
    sr1=lll$swap_rate[,1],
    sr2=lll$swap_rate[,2],
    sr3=lll$swap_rate[,3]
  )
  
  swap_rate <- rbind(swap_rate,swap_rate_temp)
  check_iter <- lll$total_iter
  if(!is_empty(check_iter)){
    iterations <-rbind(iterations,c(id_chosen,apply(check_iter, 2, mean))) 
  }
  
}


swap_rate <- swap_rate %>% 
  group_by(id) %>%
  summarise(across(everything(),mean))
saveRDS(swap_rate,"swap_rate_report.Rds")
saveRDS(iterations,"iterations_report.Rds")
print(swap_rate)

iter_report <- readRDS("results/iterations_report.Rds")
swap_rate <- readRDS("results/swap_rate_report.Rds")
View(swap_rate)
View(iter_report)
