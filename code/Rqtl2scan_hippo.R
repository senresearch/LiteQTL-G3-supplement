library(qtl)
## readin data in R/qtl
#bxd <- read.cross(file="data/hippocampus-rqtl.csv",format="csv",crosstype="risib",genotypes=c("B","D"))
print("done read.cross")

#drop obs. & traits with all NAs 
#keepidx<-which(rowSums(is.na(bxd$pheno))<1000798)

#c1<-subset(bxd,ind=keepidx)
#rownames(c1$pheno)<-c1$pheno$ID
#c1$pheno<-c1$pheno[,-1]


# ##generating a phenotype data file
# library(mice)
# library(lattice)
# #write.csv(c1$pheno,file="hippocampus-pheno.csv")
# trait<-read.csv(file="hippocampus-pheno.csv")
# end<-dim(trait)[2]
# #check NAs
# table(colSums(is.na(c1$pheno[,-end])))
# drop.idx<-which(colSums(is.na(trait[,2:(end-1)]))>42)
# trait<-trait[,2:(end-1)]
# trait<-trait[,-drop.idx]
# #write.csv(trait,file="hippocampus-pheno-nomissing.csv",row.names=FALSE)

pheno_nomissing = read.csv(file="data/hippocampus-pheno-nomissing.csv")
print("done reading pheno")

library(parallel)
library(tictoc)
library(qtl2)
#count # of cores
ncores=detectCores()
print("Number of cores:")
print(ncores)

# convert a cross from the qtl format to the qtl2 format
cvt1<-convert2cross2(c1)
#insert pseudomarker
map <- insert_pseudomarkers(cvt1$gmap, step=0)
pr <- calc_genoprob(cvt1, map, error_prob=0.002, cores=ncores)
print("done calc genoprob")

#get whole genotype prob file
getGenopr<-function(x){
  temp<<-NULL
  m=length(attributes(x)$names)
  cnames<-attributes(x)$names
  for (i in 1:m) {
    d<-eval(parse(text=paste(c('dim(x$\'', cnames[i] ,'\')'),collapse='')))
    nam<-eval(parse(text=paste(c('dimnames(x$\'',cnames[i],'\')[[2]]'),collapse = '')))
    cnam<-rep(nam,d[3])
    p_chr<-paste(c('array(x$\'',cnames[i],'\',dim=c(d[1],d[2]*d[3]))'),collapse='')
    prob<-eval(parse(text = p_chr))
    temp<-cbind(temp,prob)
  }
  return(temp)
}
prob1<-getGenopr(pr)
#write.csv(prob1,file="hippocampus-genopr-AA-BB.csv",row.names = FALSE)
print("done getting geno prob")

tic()
# out <- scan1(pr, cvt1$pheno, cores=ncores)
out <- scan1(pr, pheno_nomissing, cores=ncores)
toc()
print("done scanning, printing out result")

write.csv(out,filename="../rqtl_lod_score.csv")


