# This script is adapted from selectBXD.R
# This script cleans phenotype and genotype data of bxd family


## load library
library(tidyverse)
## read expression traits as character
d <- read.csv(file="../data/raw/bxdspleen.txt",skip=32,
              sep="\t",colClasses="character")
## read genotype file as character
g <- read.csv("../data/raw/bxd.geno",sep="\t",skip=21,
              colClasses="character")

########################
## expression traits
########################

## select probeset column
d0  <- select(d,"ProbeSet")
## select BXD columns
d1 <- select(d,matches("BXD"))
## put probeset and BXD together
d2  <- cbind(d0,d1)
## get the names of each column
d2names  <- names(d2)

########################
## genotypes
########################

## select marker information columns
g0 <- select(g,one_of("Locus","Chr","cM","Mb"))
## select BXD columns
g1 <- select(g,matches("BXD"))
## put together BXD and marker info
g2  <- cbind(g0,g1)
## get the names of each column
g2names  <- names(g2)

######################
## putting genotype and traits together
#####################

## transpose both genotype and expression data
g3  <- as_tibble(t(g2))
d3  <- as_tibble(t(d2))
## create id columns for both datasets
g3$id  <- g2names
d3$id  <- d2names

## make a right join on id
## this will keep all the traits with genotypes
gd <- right_join(d3,g3,"id")
## fill in probeset names
gd[1,1:ncol(d3)] <- d3[1,]
## create an extra ID column
ID <- as.character(gd$id)
## put that at the beginning of the data frame
gd <- cbind(ID,gd)
gd[,1] <- as.character(gd[,1])
gd[1:4,1] <- c("ID",NA,NA,NA)
## write out file
write_csv(gd,path="../data/processed/spleen-geno-pheno.csv",col_names=F,na="")

## remove mb positions
gd <- gd[-4,]
## remove id
ididx <- which(names(gd)=="id")
gd <- gd[,-ididx]

## write in R/qtl format
write_csv(gd,path="../data/processed/spleen-geno-pheno-rqtl.csv",col_names=F,na="")

## read in data in R/qtl
# bxd <- read.cross(file="../data/processed/spleen-geno-pheno-rqtl.csv",format="csv",
#                   crosstype="risib",genotypes=c("B","D"))
