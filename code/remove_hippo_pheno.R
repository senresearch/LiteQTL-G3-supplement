library(qtl)
## readin data in R/qtl
bxd <- read.cross(file="hippocampus-rqtl.csv",format="csv",crosstype="risib",genotypes=c("B","D"))

#drop obs. & traits with all NAs 
keepidx<-which(rowSums(is.na(bxd$pheno))<1000798)

c1<-subset(bxd,ind=keepidx)
rownames(c1$pheno)<-c1$pheno$ID
c1$pheno<-c1$pheno[,-1]


##generating a phenotype data file
library(mice)
library(lattice)
#write.csv(c1$pheno,file="hippocampus-pheno.csv")
# trait<-read.csv(file="hippocampus-pheno.csv")
trait <- c1$pheno
end<-dim(trait)[2]
#check NAs
table(colSums(is.na(c1$pheno[,-end])))
# TODO: double check if the 42 makes sense. 
drop.idx<-which(colSums(is.na(trait[,2:(end-1)]))>42)
trait<-trait[,2:(end-1)]
trait<-trait[,-drop.idx]
write.csv(trait,file="hippocampus-traits-nomissing.csv",row.names=FALSE)