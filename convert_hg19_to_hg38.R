#################################
#Changing Hg19 to Hg38############

# Problem: Dan's whitelist is aligned to Hg19, but our VCF files are aligned to the newer Hg38, therefore we want to convert the two
#Bioconductor has liftover which allows us to do this.
library(liftOver)
library(GenomicFeatures)
library(rtracklayer)

#Set working directory 
wd<- "/home/igw24/Documents/SGK1_data"
setwd(wd)
#########

#firstly we need to make a gRanges object containng all the SNP locations. 
whitelist<- read_csv("SGK1_whitelist.csv")
SGK1.whitelist.granges<- function(x){
                                  SGK1snps<- GRanges( seqnames= "chr6", ranges= IRanges(start = x$Pos, width =1))
                                  values(SGK1snps)<- x$Sample
                                  score(SGK1snps)<- x$"Unique ID"
                                  print(SGK1snps)}

SGK1.hg19.granges<- SGK1.whitelist.granges(whitelist)

# Get the chain file for liftOver
path =system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
ch =import.chain(path)
ch
#use liftover command
seqlevelsStyle(SGK1.hg19.granges)= "UCSC"
cur38= liftOver(SGK1.hg19.granges, ch)
class(cur38)
cur38= unlist(cur38)
genome(cur38)= "hg38"
cur38

#add the extra data
library(tidyverse)
SGK1.whitelist.hg38<- as_tibble(cur38)
SGK1.whitelist.hg38<- SGK1.whitelist.hg38[,c(7,1:6)]
colnames(SGK1.whitelist.hg38)[1]<- "UniqueID"
colnames(SGK1.whitelist.hg38)[7]<- "Sample"
colnames(whitelist)[2]<-"UniqueID"
whitelist<- as_tibble(whitelist)
SGK1.whitelist.hg38<- left_join(SGK1.whitelist.hg38, whitelist, by= c("UniqueID"))
SGK1.whitelist.hg38<- select(SGK1.whitelist.hg38, -matches("Pos"))
write.csv(SGK1.whitelist.hg38, file= "SGK1_whitelist_hg38.csv")


