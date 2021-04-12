include("new-functions-liteqtl.jl")

# Start with correlation matrix. This is just a random matrix. Nothing special


geno = [0.0 0.0 1.0 0.0; 0.0 0.0 2.0 1.0; 1.0 0.0 0.0 2.0; 0.0 1.0 0.0 0.0]

# pheno = [1.9 1.2 9.3 2.4; 8.5 3.6 9.7 2.8; 0.6 2.4 3.3 4.5; 0.1 6.7 5.8 9.2]

G = filter_maf(geno)


# # Degree of freedom should be n - 1
# mydof = size(geno, 1) -1
# n = size(pheno,1)
# pheno_std = LiteQTL.get_standardized_matrix(pheno)
# geno_std = LiteQTL.get_standardized_matrix(geno)

# r = LiteQTL.calculate_r(pheno_std, geno_std)
# lodc = LiteQTL.lod_score_multithread(n, r )
# pval = lod2p(lodc)

# r = calculate_r_tensor(geno, pheno)
# pval = pval_calc(r, mydof)
# display(pval)
