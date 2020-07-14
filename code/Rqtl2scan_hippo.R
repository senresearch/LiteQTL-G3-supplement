# This script is adapted from Rqtl2scan_hippo.R

library(qtl)
## readin data in R/qtl
bxd <- read.cross(file="../data/processed/hippocampus-geno-pheno-rqtl.csv",format="csv",crosstype="risib",genotypes=c("B","D"))
print("done read.cross")

#drop obs. & traits with all NAs
keepidx<-which(rowSums(is.na(bxd$pheno))<1000798)

c1<-subset(bxd,ind=keepidx)
rownames(c1$pheno)<-c1$pheno$ID
c1$pheno<-c1$pheno[,-1]


trait<-c1$pheno
end<-dim(trait)[2]

 #check NAs
table(colSums(is.na(trait[,-end])))
drop.idx<-which(colSums(is.na(trait[,2:(end-1)]))>42)
trait<-trait[,1:(end-1)]
trait<-trait[,-drop.idx]
write.csv(trait, file="../data/processed/hippo-pheno-nomissing.csv")


library(parallel)
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
write.csv(pr, file="../data/processed/bxd-genoprob_hippo.csv")

scantime <- system.time({
    out <- scan1(pr, trait, cores=ncores)
  })
print("Rqtl genome scan for hippocampus data took ", scantime)
print("done scanning, printing out result...")

write.csv(out,file="../data/results/rqtl_lod_score_hippo.csv")
