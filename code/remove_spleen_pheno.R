library(qtl)
## readin data in R/qtl
bxd <- read.cross(file="../data/processed/spleen-geno-pheno-rqtl.csv",format="csv",
                  crosstype="risib",genotypes=c("B","D"))



#drop obs. & traits with all NAs 
keepidx<-which(rowSums(is.na(bxd$pheno))<35500)

c1<-subset(bxd,ind=keepidx)
rownames(c1$pheno)<-c1$pheno$ID
c1$pheno<-c1$pheno[,-1]

# droptrait<-which(colSums(is.na(c1$pheno)) == 79) # From Hyeonju's code. But it didn't remove all NAs. 
droptrait<-which(colSums(is.na(c1$pheno))>0)
c1$pheno<-c1$pheno[,-droptrait]


end<-dim(c1$pheno)[2]
c1$pheno <- c1$pheno[,1:end-1]



write.csv(c1$pheno, file="../data/processed/spleen_traits_nomissing.csv")
