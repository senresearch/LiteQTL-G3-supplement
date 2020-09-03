
## load library
library(tidyverse)
library(data.table)

## read expression traits as character

# > dim(pheno1)
# [1] 1000798     108
# > dim(pheno2)
# [1] 356054    108
# Warning message:
# In fread(file = "../data/raw/bxdhippo.txt", skip = 32, sep = "\t",  :
#   Stopped early on line 356089. Expected 108 fields but found 9. Consider fill=TRUE and comment.char=. First discarded non-empty line: <<4597261      Phr1    14      103.51306     105689  -       103.51306       103.513293      pam, highwire, rpm 1>>

reorgdata<-function(genofile, phenofile, outputfile){

  pheno <- read.csv(phenofile,skip=32,
                sep="\t",colClasses="character")

  ## read genotype file as character
  geno <- read.csv(genofile,sep="\t",skip=21,
                colClasses="character")


  ########################
  ## expression traits
  ########################

  ## select probeset column
  probeset  <- select(pheno,"ProbeSet")
  ## select BXD columns
  pheno_bxd_cols <- select(pheno,matches("BXD"))
  ## put probeset and BXD together
  chosenpheno  <- cbind(probeset,pheno_bxd_cols)
  ## get the names of each column
  chosenphenonames  <- names(chosenpheno)

  ########################
  ## genotypes
  ########################

  ## select marker information columns
  gmap <- select(geno,one_of("Locus","Chr","cM","Mb"))
  ## select BXD columns
  geno_bxd_cols <- select(geno,matches("BXD"))
  ## put together BXD and marker info
  chosengeno  <- cbind(gmap,geno_bxd_cols)
  ## get the names of each column
  chosengenonames  <- names(chosengeno)

  ######################
  ## putting genotype and traits together
  #####################

  ## transpose both genotype and expression data
  ## make the first column the variable names of the tibble
  chosengeno_trans  <- t(chosengeno[,-1])
  colnames(chosengeno_trans) <- chosengeno[,1]
  chosengeno_tb <- as_tibble(chosengeno_trans)
  #
  chosenpheno_trans  <- t(chosenpheno[,-1])
  colnames(chosenpheno_trans)  <- chosenpheno[,1]
  chosenpheno_tb  <- as_tibble(chosenpheno_trans) 

  ## create id columns for both datasets
  chosengeno_tb$id  <- chosengenonames[-1]
  chosenpheno_tb$id  <- chosenphenonames[-1]

  ## make a right join on id
  ## this will keep all the traits with genotypes
  genopheno <- right_join(chosenpheno_tb,chosengeno_tb,"id")

  ## get the rows with the marker info; they don't start with "B"
  genophenoMkinfo <- filter(genopheno,!str_detect(id,"^B+."))
  ## all other rows
  genophenoinfo <- filter(genopheno,str_detect(id,"^B+."))
  ## bind rows and make tibble
  genopheno <- tibble(bind_rows(genophenoMkinfo,genophenoinfo))
  ## sanitize the id column by getting rid of the marker info annotation
  genopheno$id[!str_detect(genopheno$id,"^B+.")] <- ""
  ## make id the first column 
  genopheno <- relocate(genopheno,id,.before=1)

  ## remove mb positions
  genopheno <- filter(genopheno,id!="Mb")

  transposedf<-function(df){
    library(data.table)
  #   rownames(df) <- df[,1]
    t_df <- transpose(df)
  #   colnames(t_test) <- t_test[1, ]
  #   t_test <- t_test[-1, ]
    rownames(t_df) <- colnames(df)
    colnames(t_df) <- t_df[1, ] 

    return(t_df[-1, ])
  }

  t_genopheno = transposedf(genopheno)
  ## write in R/qtl format
  write_csv(t_genopheno,path=outputfile,col_names=F, na="")
}

# running genome scan for spleen data. 
genofile <- "../data/raw/bxd.geno" 
phenofile <- "../data/raw/bxdspleen.txt"
outputfile <- "../data/processed/spleen-geno-pheno-rqtl.csv"

reorgdata(genofile, phenofile, outputfile)


# running genome scan for hippocampus data. 
# genofile <- "../data/raw/bxd.geno"
# phenofile <- "../data/raw/bxdhippo.txt"
# outputfile <- "../data/processed/hippo-geno-pheno-rqtl.csv"

# reorgdata(genofile, phenofile, outputfile)

