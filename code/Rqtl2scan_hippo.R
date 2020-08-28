# This script is adapted from Rqtl2scan_hippo.R

library(qtl)
## readin data in R/qtl
bxd <- read.cross(file="../data/processed/hippo-geno-pheno-rqtl.csv",format="csvr",
                  crosstype="risib",genotypes=c("B","D"))
print("done read.cross. ")

#drop obs. & traits with all NAs
keepidx<-which(rowSums(is.na(bxd$pheno))<1000798)

c1<-subset(bxd,ind=keepidx)
end<-dim(c1$pheno)[2]

 #check NAs
table(colSums(is.na(c1$pheno)))
drop.idx<-which(colSums(is.na(c1$pheno))>0)
c1$pheno<-c1$pheno[,-drop.idx]


write.csv(t(c1$pheno), file="../data/processed/hippo-pheno-nomissing.csv", row.names=FALSE)


library(parallel)
library(qtl2)

# convert a cross from the qtl format to the qtl2 format
cvt1<-convert2cross2(c1)
#insert pseudomarker
map <- insert_pseudomarkers(cvt1$gmap, step=0)
prtime <- system.time({
    pr <- calc_genoprob(cvt1, map, error_prob=0.002, cores=4)
})
# > scantime2
#     user   system  elapsed
#   67.856 6154.076 2830.997

print("done calc genoprob")
write.csv(pr, file="../data/processed/hippo-bxd-genoprob.csv", row.names=FALSE)
# if read in, length(pr) = 14642

scantime <- system.time({
    out <- scan1(pr, cvt1$pheno, cores=1)
})
# > scantime
#     user   system  elapsed
# 12496.70 12643.84 25154.36

print("Rqtl genome scan for hippocampus data took ", scantime)
print("done scanning, printing out result...")

library(Rfast)

rowmaxidx = rowMaxs(out, value = FALSE)
rowmaxval = rowMaxs(out, value = TRUE)
write.csv(c(rowmaxidx, rowmaxval), file="../data/results/hippo_rqtl_lod_score.csv", row.names=FALSE)