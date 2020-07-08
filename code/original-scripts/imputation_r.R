#calculate genotype probability (filling genotype prob only)
p_gprob<-calc.genoprob(bxd)
#pull genotype and phenotype data
gprob<-pull.genoprob(p_gprob)
pheno<-pull.pheno(bxd)
#get indices of rows having too many NA
omitind<-which(rowSums(is.na(pheno))<35500)

pheno0<-pheno[omitind,]
geno0<-gprob[omitind,]

write.csv(geno0,file="geno_prob.csv",quote = FALSE)

pheno1<-pheno0[,-1]
rownames(pheno1)<-pheno0[,1]
#check which columns have missing values
table(colSums(is.na(pheno1)))

#get indices of colmns being empty
idx79<-which(colSums(is.na(pheno1))==79)
pheno2<-pheno1[,-idx79]

#pick random indices for imputation
pickidx<-sample(setdiff(1:dim(pheno2)[2],idx2),30)
ordix<-sort(c(pickidx,idx2))

#imputation
library(mice)
#library(lattice)

imp<-pheno2[,ordix]

imp1<-mice(imp,defaultMethod = "pmm",seed = 10)
imp2<-complete(imp1)
pheno2[,idx2[1]]<-imp2[,which(names(imp2)==names(idx2)[1])]
pheno2[,idx2[2]]<-imp2[,which(names(imp2)==names(idx2)[2])]

write.csv(pheno2,file = "traits.csv",quote = FALSE)
# write.csv(pheno3,file = "traits_with_missing.csv",quote = FALSE) : for reference data