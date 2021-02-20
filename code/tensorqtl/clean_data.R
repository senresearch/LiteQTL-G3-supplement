library(stringr)
library(parallel)
library(data.table)
library(rres)


datacleaning <- function(genotypefile, phenotypefile){
    print("Reading in genotype...")
    geno = fread(genotypefile)
    # > geno[1:5, 1:5]
    #                    snp HG00096 HG00097 HG00099 HG00100
    # 1: chr1_10177_A_AC_b38       1       1       1       1
    # 2: chr1_10352_T_TA_b38       1       1       1       1
    # 3:  chr1_11008_C_G_b38       0       0       0       0
    # 4:  chr1_11012_C_G_b38       0       0       0       0
    # 5:  chr1_13110_G_A_b38       0       1       0       0

    print("Reading in phenotype data...")
    pheno = fread(phenotypefile)
    # > pheno[1:5, 1:5]
    #    #chr  start    end            gene_id    HG00096
    # 1: chr1  29552  29553  ENSG00000227232.5 -0.3859259
    # 2: chr1 135894 135895  ENSG00000268903.1  0.7357491
    # 3: chr1 137964 137965  ENSG00000269981.1  0.5295665
    # 4: chr1 195410 195411  ENSG00000279457.4  0.3379166
    # 5: chr1 522927 522928 ENSG00000237094.12  0.5621835

    ########################
    ## expression traits
    ########################
    # no need to find missing data since there is no missing data
    # no need to match individual in pheno and geno since all individual exists in both files. 
    # removing the first three columns in phenotype, only keeping gene_id and indiv data
    pheno = pheno[, 4:dim(pheno)[2]]
    # renaming the name of the first column:
    colnames(pheno)[1] = "ID"
    fwrite(pheno, file="cleanpheno.csv", row.name=FALSE)

    ########################
    ## genotypes
    ########################
    # split the snp string to get chromosome number. 
    elems <- unlist( strsplit(geno$snp , "_" ) )
    # turn it into a matrix, we know there are 5 elements in the split string. the first element is the chromsome number. 
    m <- matrix( elems , ncol = 5 , byrow = TRUE )
    chr = m[,1]
    # recombine the genotype. 
    geno <- cbind(geno$snp, chr, geno[,2:dim(geno)[2]])
    colnames(geno)[1:2] = c("ID", "")
    # df[df$aged <= df$laclen, ] 
    chosenchr = "chr9"
    subgeno <- geno[geno$snp == chosenchr, ]
    fwrite(subgeno, file="chr9.csv", row.names=FALSE)
}

genotypefile = "../../data/tensorqtldata/genotype.ped"
phenotypefile = "../../data/tensorqtldata/GEUVADIS.445_samples.expression.bed"

datacleaning(genotypefile, phenotypefile)





