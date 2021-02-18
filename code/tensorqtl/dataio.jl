using CSV 
# reading genotype:
geno = CSV.read("genotype.csv")
# reading phenotype: # note to self: need to gunzip the BED file first. 
pheno = CSV.read("GEUVADIS.445_samples.expression.bed")
# read covariates:
covar = CSV.read("GEUVADIS.445_samples.covariates.txt")

pheno[1:5, 1:10]
geno[1:5, 1:10]
covar[1:5, :]
println("$(size(geno)), $(size(pheno)), $(size(covar))")

using DelimitedFiles 
bxdgeno = readdlm("/home/xiaoqihu/git/lmgpu_bin/data/raw/bxd.geno", skipstart=21)

