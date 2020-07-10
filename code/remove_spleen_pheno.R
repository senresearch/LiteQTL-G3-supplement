library(qtl)
## readin data in R/qtl
bxd <- read.cross(file="../data/processed/spleen-geno-pheno-rqtl.csv",format="csv",
                  crosstype="risib",genotypes=c("B","D"))



#drop obs. & traits with all NAs 
keepidx<-which(rowSums(is.na(bxd$pheno))<35500)

c1<-subset(bxd,ind=keepidx)
rownames(c1$pheno)<-c1$pheno$ID
c1$pheno<-c1$pheno[,-1]

droptrait<-which(colSums(is.na(c1$pheno))==79)
c1$pheno<-c1$pheno[,-droptrait]

write.csv(c1$pheno, file="../data/processed/spleen_traits_nomissing.csv")
# pheno<-read.csv("traits.csv",sep=",")
# c1$pheno<-pheno
