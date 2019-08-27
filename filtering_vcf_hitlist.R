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

# using the dplyr function

sgk1hitlist.tbl<- as_tibble(SGK1_hitlist)


xy<-sgk1hitlist.tbl %>% filter(FILTER == "PASS") %>%
                        mutate(confidence= case_when(
                          depthALT<8 ~"low confidence",
                          CONSEQUENCE== "nonsynonymous" &  depthALT>8 & start%in%SGK1_whitelist_hg38$start ~"High Confidence",
                          CONSEQUENCE== "synonymous" &depthALT>8 & start%in%SGK1_whitelist_hg38$start   ~"medium confidence" ,
                          depthALT>8 & !start%in%SGK1_whitelist_hg38$start ~ "not in whitelist",
                          TRUE ~"other"))

xy%>% count(confidence)

     