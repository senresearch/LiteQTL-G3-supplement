# This script is adapted from code/original_scripts/Rqtl2scan.R

library(qtl)
## readin data in R/qtl
# bxd <- read.cross(file="../data/processed/spleen-geno-pheno-rqtl.csv",format="csv",
#                   crosstype="risib",genotypes=c("B","D"))

bxd <- read.cross(file="../data/processed/spleen-geno-pheno-rqtl.csv",format="csv",
                  crosstype="risib",genotypes=c("B","D"))

#  --Read the following data:
#          198  individuals
#          7321  markers
#          35559  phenotypes
#  --Cross type: risib
# Warning messages:
# 1: In read.cross.csv(dir, file, na.strings, genotypes, estimate.map,  :
#   The following unexpected genotype codes were treated as missing.
#     |H|

# 2: In summary.cross(cross) :
#   Duplicate markers [1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 2, 3, 4, 5, 6, 7, 8, 9, X]

#drop obs. & traits with all NAs 
keepidx<-which(rowSums(is.na(bxd$pheno))<35500)

c1<-subset(bxd,ind=keepidx)
rownames(c1$pheno)<-c1$pheno$ID
c1$pheno<-c1$pheno[,-1]


end<-dim(c1$pheno)[2]


droptrait<-which(colSums(is.na(c1$pheno))>0)
c1$pheno<-c1$pheno[,1:(end-1)]
c1$pheno<-c1$pheno[,-droptrait]


write.csv(t(c1$pheno), file="../data/processed/spleen-pheno-nomissing.csv", row.names=FALSE)
# extract genotype data from the processed data
#gen<-pull.geno(c1)
#write.csv(gen,file="genotypedata.csv")
# pheno<-read.csv("traits.csv",sep=",")
# c1$pheno<-pheno

library(parallel)
library(qtl2)

# convert a cross from the qtl format to the qtl2 format
cvt1<-convert2cross2(c1)


#insert pseudomarker
map <- insert_pseudomarkers(cvt1$gmap, step=0)
prtime <- system.time({
    pr <- calc_genoprob(cvt1, map, error_prob=0.002, cores=4)
})

write.csv(pr, file="../data/processed/bxd-genoprob_spleen.csv", row.names=FALSE)

scantime <- system.time({
    out <- scan1(pr, cvt1$pheno, cores=4)
})
# > scantime 
#    user  system elapsed
#  29.956 445.484 180.510

# > scantime # using cores=1
#    user  system elapsed
# 461.284  10.776 472.084

cat("Rqtl genome scan for spleen data took ", scantime)
# This ran 56.461 seconds on imac pro (3 GHz Intel Xeon W)
write.csv(out,file="../data/results/rqtl_lod_score_spleen.csv",row.names=FALSE)

