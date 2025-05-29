#!/usr/bin/env Rscript
library("optparse")
option_list = list(
  make_option(c("-i", "--init"), type="integer", default=1, 
              help="Initial value set [default= %default]", metavar="integer"),
  make_option(c("-c", "--constant"), type="numeric", default=2, 
              help="geometric constant [default= %default]", metavar="numeric"),
  make_option(c("-u", "--upper_index"), type="numeric", default=10, 
              help="Upper bound for geometric schedule according to highest power [default= %default]", metavar="numeric"),
  make_option(c("-t", "--ntemp"), type="integer", default=10, 
              help="Number of temperatures [default= %default]", metavar="integer"),
  make_option(c("-n", "--nsim"), type="integer", default=10000, 
              help="Number of samples [default= %default]", metavar="integer"),
  make_option(c("-l", "--logit"), type="character", default="y", 
              help="Transform transition probabilities using multinomial logit function [default= %default]", metavar="character"),
  make_option(c("-a", "--auxiliar"), type="character", default="n", 
              help="Use latent variables for transition probabilities [default= %default]", metavar="character"),
  make_option(c("-x", "--within"), type="integer", default="1", 
              help="Number of iterations before swap proposal [default= %default]", metavar="integer"),
  make_option(c("-b", "--batch"), type="integer", default="1", 
              help="1 million batch number [default= %default]", metavar="integer"),
  make_option(c("-d", "--sdt"), type="character", default="n", 
              help="Step size tunning dependent on temperature [default= %default]", metavar="character"),
  make_option(c("-s", "--stimulli"), type="character", default="n", 
              help="Include sound stimulli in transition probabilities [default= %default]", metavar="character"),
  make_option(c("-e", "--emission"), type="integer", default=0, 
              help="Emission stream about to explored [default= %default]", metavar="integer"),
  make_option(c("-v", "--version"), type="integer", default=1, 
              help="Temperature schedule version [default= %default]", metavar="integer")
  
  
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

sourceCpp("you_Cpp_file")

# Here you extract the flags provided when running Rscript
temp = opt$ntemp
nsimm = opt$nsim
opt$log_emissions = "y"
explore_in_log = opt$log_emissions == "y"
logit = opt$logit == "y"
latent_var = opt$auxiliar == "y"
sdt = opt$sdt == "y"
within_temp = opt$within
include_covariate = opt$stimulli == "y"
