# This script is adapted from Rqtl2scan_hippo.R

library(qtl)
## readin data in R/qtl
bxd <- read.cross(file="../data/processed/hippo-geno-pheno-rqtl.csv",format="csvr",
                  crosstype="risib",genotypes=c("B","D"))
print("done read.cross. ")

#drop obs. & traits with all NAs
keepidx<-which(rowSums(is.na(bxd$pheno))<1000798)

c1<-subset(bxd,ind=keepidx)
# rownames(c1$pheno)<-c1$pheno$ID
# c1$pheno<-c1$pheno[,-1]


trait<-c1$pheno
end<-dim(trait)[2]

 #check NAs
# trait<-trait[,1:(end-1)]
table(colSums(is.na(trait)))
drop.idx<-which(colSums(is.na(trait))>0)
trait<-trait[,-drop.idx]

# transposedf<-function(df){
#   library(data.table)
#   t_df <- transpose(df)
#   colnames(t_df) <- rownames(df)
#   rownames(t_df) <- colnames(df)
#   return(t_df)
# }

# transposedf<-function(df){
#   library(data.table)
#   t_df <- transpose(df)
#   rownames(t_df) <- colnames(df)
#   colnames(t_df) <- t_df[1, ] 
#   return(t_df[-1, ])
# }


write.csv(t(trait), file="../data/processed/hippo-pheno-nomissing.csv", row.names=FALSE)


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

# print("done calc genoprob")
write.csv(pr, file="../data/processed/bxd-genoprob_hippo.csv", row.names=FALSE)

scantime <- system.time({
    out <- scan1(pr, cvt1$pheno, cores=1)
})
# > scantime
#     user   system  elapsed
# 12496.70 12643.84 25154.36

print("Rqtl genome scan for hippocampus data took ", scantime)
print("done scanning, printing out result...")

write.csv(out,file="../data/results/rqtl_lod_score_hippo.csv", row.names=FALSE)