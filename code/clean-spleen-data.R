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
## make the first column the variable names of the tibble
g3  <- t(g2[,-1])
colnames(g3) <- g2[,1]
g4 <- as_tibble(g3)
#
d3  <- t(d2[,-1])
colnames(d3)  <- d2[,1]
d4  <- as_tibble(d3)

## create id columns for both datasets
g4$id  <- g2names[-1]
d4$id  <- d2names[-1]

## make a right join on id
## this will keep all the traits with genotypes
gd <- right_join(d4,g4,"id")

## get the rows with the marker info; they don't start with "B"
gdMkinfo <- filter(gd,!str_detect(id,"^B+."))
## all other rows
gdinfo <- filter(gd,str_detect(id,"^B+."))
## bind rows and make tibble
gd <- tibble(bind_rows(gdMkinfo,gdinfo))
## sanitize the id column by getting rid of the marker info annotation
gd$id[!str_detect(gd$id,"^B+.")] <- ""
## make id the first column 
gd <- relocate(gd,id,.before=1)

## write out file
write_csv(gd,path="../data/processed/spleen-geno-pheno.csv",col_names=F,na="")

## remove mb positions
gd <- filter(gd,id!="Mb")

## write in R/qtl format
write_csv(gd,path="../data/processed/spleen-geno-pheno-rqtl.csv",col_names=F,na="")

## read in data in R/qtl
# bxd <- read.cross(file="../data/processed/spleen-geno-pheno-rqtl.csv",format="csv",
#                  crosstype="risib",genotypes=c("B","D"))
