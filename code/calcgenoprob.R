# library(qtl)
# ## readin data in R/qtl
# bxd <- read.cross(file="../data/processed/spleen-geno-pheno-rqtl.csv",format="csv",crosstype="risib",genotypes=c("B","D"))
# 
# #drop obs. & traits with all NAs 
# keepidx<-which(rowSums(is.na(bxd$pheno))<35500)
# 
# c1<-subset(bxd,ind=keepidx)
# library(parallel)
# 
# library(qtl2)
# #count # of cores
# ncores=detectCores()
# # convert a cross from the qtl format to the qtl2 format
# cvt1<-convert2cross2(c1)
# #insert pseudomarker
# map <- insert_pseudomarkers(cvt1$gmap, step=1)
# pr <- calc_genoprob(cvt1, map, error_prob=0.002, cores=ncores)
# 
# #get whole genotype prob file
# getGenopr<-function(x){
#   temp<<-NULL
#   m=length(attributes(x)$names)
#   cnames<-attributes(x)$names
#   for (i in 1:m) {
#     d<-eval(parse(text=paste(c('dim(x$\'', cnames[i] ,'\')'),collapse='')))
#     nam<-eval(parse(text=paste(c('dimnames(x$\'',cnames[i],'\')[[2]]'),collapse = '')))
#     cnam<-rep(nam,d[3])
#     p_chr<-paste(c('array(x$\'',cnames[i],'\',dim=c(d[1],d[2]*d[3]))'),collapse='')
#     prob<-eval(parse(text = p_chr))
#     temp<-cbind(temp,prob)
#   }
#   return(temp)
# }
# prob1<-getGenopr(pr)
# write.csv(prob1,file="../data/processed/genopr-AA-BB.csv",row.names = FALSE)


# ######### Got this script from imputation_r.R
library(qtl)
# ## readin data in R/qtl

bxd <- read.cross(file="../data/processed/spleen-geno-pheno-rqtl.csv",format="csv",
                  crosstype="risib",genotypes=c("B","D"))
# calculate genotype probability (filling genotype prob only)
p_gprob<-calc.genoprob(bxd)
#pull genotype and phenotype data
gprob<-pull.genoprob(p_gprob)
pheno<-pull.pheno(bxd)

# #get indices of rows having too many NA
omitind<-which(rowSums(is.na(pheno))<35500)


geno0<-gprob[omitind,]

write.csv(geno0,file="../data/processed/geno_prob.csv",quote = FALSE)

