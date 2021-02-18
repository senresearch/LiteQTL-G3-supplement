library(qtl)
library(stringr)
library(parallel)
library(qtl2)
library(data.table)

timefunc <- function(message, func, ...){
  t = system.time({
    val = func(...)
  })
  print("--------")
  print(t)
  print("--------")
  return (val)
}
getmatchingBXD <- function(phenonames, genonames){
  # find all the matching columnnames in pheno and geno. 
  matched = phenonames[which(phenonames %in% genonames)]
  # filter out BXD columns. They start with B
  return (matched[str_detect(matched,"^B+.")])
}

datacleaning <- function(genotypefile, phenotype){
    print("Reading in genotype...")
    geno = fread(genotypefile)
    # > geno[1:5, 1:5]
    #                    snp HG00096 HG00097 HG00099 HG00100
    # 1: chr1_10177_A_AC_b38       1       1       1       1
    # 2: chr1_10352_T_TA_b38       1       1       1       1
    # 3:  chr1_11008_C_G_b38       0       0       0       0
    # 4:  chr1_11012_C_G_b38       0       0       0       0
    # 5:  chr1_13110_G_A_b38       0       1       0       0

    print("Reading in phenotype data...")
    pheno = fread(phenotypefile)
    # > pheno[1:5, 1:5]
    #    #chr  start    end            gene_id    HG00096
    # 1: chr1  29552  29553  ENSG00000227232.5 -0.3859259
    # 2: chr1 135894 135895  ENSG00000268903.1  0.7357491
    # 3: chr1 137964 137965  ENSG00000269981.1  0.5295665
    # 4: chr1 195410 195411  ENSG00000279457.4  0.3379166
    # 5: chr1 522927 522928 ENSG00000237094.12  0.5621835

    ########################
    ## expression traits
    ########################
    # no need to find missing data since there is no missing data
    # no need to match individual in pheno and geno since all individual exists in both files. 
    # removing the first three columns in phenotype, only keeping gene_id and indiv data
    pheno = pheno[, 4:dim(pheno)[2]]
    # renaming the name of the first column:
    colnames(pheno)[1] = "ID"
    fwrite(pheno, file="cleanpheno.csv", row.name=FALSE)

    ########################
    ## genotypes
    ########################
    chr = sapply(geno[,1], substr, 4,4)
    geno <- cbind(geno$snp, chr, geno[,2:dim(geno)[2]])
    colnames(geno)[1:2] = c("ID", "")
    fwrite(geno, file="cleangeno.csv", row.names=FALSE)

    bxd <- read.cross("csvsr", ".", rqtlgenofile, cleanphenofile, crosstype = "risib", genotypes = c("B", "D"))

}
reorgdata<-function(genofile, phenofile, cleanphenofile, genoprobfile, gmapfile){
  
  pheno <- read.csv(phenofile,skip=32,sep="\t",colClasses="character",na.strings="")
  print("Done reading phenotype. ")
  ## read genotype file as character
  geno <- read.csv(genofile,sep="\t",skip=21,colClasses="character",na.strings="")
  print("Done reading genotype. ")
  rqtlgenofile="temprqtlgeno.csv"
  # rqtlphenofile="./data/processed/temprqtlpheno.csv"

  ########################
  ## expression traits
  ########################
  # remove individuals that has all NAs. 
  emptyindiv <- which(colSums(is.na(pheno))>nrow(pheno)-1)
  pheno <- pheno[,-emptyindiv]
  # find individuals that has data in both pheno and geno. 
  matchingbxds <- getmatchingBXD(names(pheno), names(geno))
  phenobxd <- pheno[,matchingbxds]
  # get probeset column, this will be used as IDs. Then combine with BXD columbs. 
  subpheno <- cbind(pheno$ProbeSet, phenobxd)
  #change "ProbeSet" to "ID", this is to match rqtl format, ID column must match with genotype
  colnames(subpheno)[1] = "ID"
  # remove individuals that has missing data. 
  drop.idx<-which(rowSums(is.na(subpheno))>0)
  subpheno<-subpheno[-drop.idx,]
  write.csv(subpheno, file=cleanphenofile, row.names=FALSE)

  ########################
  ## genotypes
  ########################z
  genobxd <- geno[,matchingbxds]
  gmap <- geno[,c("Locus","Chr","cM","Mb")]
  write.csv(gmap, file=gmapfile, na="")
  # get probeset column, this will be used as IDs. 
  subgeno <- cbind(geno$Locus, geno$Chr, geno$cM, geno[,matchingbxds])
  # changed columns
  colnames(subgeno)[1:3] = c("ID", "", "")
  write.csv(subgeno, file=rqtlgenofile, row.names=FALSE)

  bxd <- read.cross("csvsr", ".", rqtlgenofile, cleanphenofile, crosstype = "risib", genotypes = c("B", "D"))

  if (file.exists(rqtlgenofile)) {file.remove(rqtlgenofile)}
  # if (file.exists(rqtlphenofile)) {file.remove(rqtlphenofile)}

  cvt1<-convert2cross2(bxd)
  print(paste("Number of chromosomes is", n_chr(cvt1)))
  map <- insert_pseudomarkers(cvt1$gmap, step=0)
  pr <- calc_genoprob(cvt1, map, error_prob=0.002, cores=1)
  print("done calc genoprob")
  write.csv(pr, file=genoprobfile, row.names=FALSE)
}

args = commandArgs(trailingOnly=TRUE)
reorgdata(args[1], args[2], args[3], args[4], args[5])


