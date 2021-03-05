scan <- function(rqtlfile, cleanphenofile, genoprobfile, rqtlresultfile, maxlod){
    library(qtl)
    bxd <- read.cross(file=rqtlfile,format="csvr",
                    crosstype="risib",genotypes=c("B","D"))
    print(paste("Number of chromosomes is", nchr(bxd)))
    print("done read.cross. ")
    #drop obs. & traits with all NAs
    keepidx<-which(rowSums(is.na(bxd$pheno))<ncol(bxd$pheno)-1)
    c1<-subset(bxd,ind=keepidx)
    end<-dim(c1$pheno)[2]

    #check NAs
    table(colSums(is.na(c1$pheno)))
    drop.idx<-which(colSums(is.na(c1$pheno))>0)
    c1$pheno<-c1$pheno[,-drop.idx]
    write.csv(c1$pheno, file=cleanphenofile, row.names=TRUE)

    library(parallel)
    library(qtl2)
    cvt1<-convert2cross2(c1)
    print(paste("Number of chromosomes is", n_chr(cvt1)))
    map <- insert_pseudomarkers(cvt1$gmap, step=0)
    prtime <- system.time({
        pr <- calc_genoprob(cvt1, map, error_prob=0.002, cores=1)
    })
    print("Calculating genoprob took: ")
    print(prtime)

    print("done calc genoprob")
    write.csv(pr, file=genoprobfile, row.names=FALSE)
    write.csv(cvt1$gmap, file=gmapfile, row.names=FALSE)

    scantime <- system.time({
        out <- scan1(pr, cvt1$pheno, cores=16)
    })

    print("Rqtl genome scan took: ")
    print(scantime)

    if (maxlod) {
        colmaxval = c()
        colmaxidx = c()
        for (i in 1:ncol(out)) {
            # colmaxval[i] = max(out[, i])
            colmaxidx[i] = which.max(out[, i])
            colmaxval[i] = out[colmaxidx[i], i]
        }
        out <- cbind(colmaxidx, colmaxval)
        print("writing max lod...")
    }
    write.csv(out, file=rqtlresultfile, row.names=FALSE)

}

# Run scan for spleen data. 
gmapfile <- "../data/processed/gmap.csv"

rqtlfile <- "../data/processed/spleen-geno-pheno-rqtl.csv"
cleanphenofile <- "../data/processed/spleen-pheno-nomissing.csv"
genoprobfile <- "../data/processed/spleen-bxd-genoprob.csv"
rqtlresultfile <- "../data/results/spleen_rqtl_lod_score.csv"
maxlod = TRUE

scan(rqtlfile, cleanphenofile, genoprobfile, rqtlresultfile, maxlod)

## Run scan for hippo data. 

rqtlfile <- "../data/processed/hippo-geno-pheno-rqtl.csv"
cleanphenofile <- "../data/processed/hippo-pheno-nomissing.csv"
genoprobfile <- "../data/processed/hippo-bxd-genoprob.csv"
rqtlresultfile <- "../data/results/hippo_rqtl_lod_score.csv"
maxlod = TRUE 

scan(rqtlfile, cleanphenofile, genoprobfile, rqtlresultfile, maxlod)
