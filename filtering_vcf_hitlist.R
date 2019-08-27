#######################################
#####Filtered the VCF hitlist file#####

library(tidyverse)

#Set working directory 
wd<- "/home/igw24/Documents/SGK1_data"
setwd(wd)
#########

#get the files
SGK1_whitelist_hg38 <- read_csv("SGK1_whitelist_hg38.csv")
SGK1_hitlist <- read_csv("SGK1_mutations2.csv")
