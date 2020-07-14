# This script is adapted from code/original_scripts/Rqtl2scan.R

library(qtl)
## readin data in R/qtl
bxd <- read.cross(file="../data/processed/spleen-geno-pheno-rqtl.csv",format="csv",
                  crosstype="risib",genotypes=c("B","D"))



#drop obs. & traits with all NAs 
keepidx<-which(rowSums(is.na(bxd$pheno))<35500)

c1<-subset(bxd,ind=keepidx)
rownames(c1$pheno)<-c1$pheno$ID
c1$pheno<-c1$pheno[,-1]


end<-dim(c1$pheno)[2]


droptrait<-which(colSums(is.na(c1$pheno))>0)
c1$pheno<-c1$pheno[,1:(end-1)]
c1$pheno<-c1$pheno[,-droptrait]
write.csv(c1$pheno, file="../data/processed/spleen-pheno-nomissing.csv")
# extract genotype data from the processed data
#gen<-pull.geno(c1)
#write.csv(gen,file="genotypedata.csv")
# pheno<-read.csv("traits.csv",sep=",")
# c1$pheno<-pheno

library(qtl2)
# convert a cross from the qtl format to the qtl2 format
cvt1<-convert2cross2(c1)
#insert pseudomarker
map <- insert_pseudomarkers(cvt1$gmap, step=0)
pr <- calc_genoprob(cvt1, map, error_prob=0.002, cores=4)
write.csv(pr, file="../data/processed/bxd-genoprob_spleen.csv")

scantime <- system.time({
    out <- scan1(pr, cvt1$pheno, cores=32)
})
cat("Rqtl genome scan for spleen data took ", scantime)
# This ran 56.461 seconds on imac pro (3 GHz Intel Xeon W)
write.csv(out,file="../data/results/rqtl_lod_score_spleen.csv")

