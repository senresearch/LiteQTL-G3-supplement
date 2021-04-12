import numpy as np
import torch 
from scipy import stats
import statistics 


def pval_calc(r, dof):
    tstat_t = r * torch.sqrt(dof / (1 - r**2))
    # tstatcpu = tstat_t.cpu()
    pval = 2*stats.t.cdf(-np.abs(tstat_t.cpu()), dof)
    return pval 

def center_normalize(matrix, dim=1):
    newmatrix = matrix - matrix.mean(dim=dim, keepdim=True)
    return newmatrix /torch.sqrt(torch.pow(newmatrix, 2).sum(dim=dim, keepdim=True))

def calculate_corr(genotype, phenotype):
    genotype = center_normalize(genotype, dim=1)
    phenotype = center_normalize(phenotype, dim=1)
    return torch.mm(genotype, phenotype.t())

def p_corr(df1, df2):
    corr = df1.corr(df2)
    N = np.sum(df1.notnull())
    t = corr*np.sqrt((N-2)/(1-corr**2))
    p = 1-scipy.stats.t.cdf(abs(t),N-2)  # one-tailed
    return corr, t, p 

def filter_maf(genotype_t, variant_id, maf_threshold):
    # """Calculate minor allele frequency"""
    alleles = 2
    af_t = genotype_t.sum(1) / (alleles * genotype_t.shape[1])
    maf_t = torch.where(af_t > 0.5, 1 - af_t, af_t)
    # """Calculate MAF and filter genotypes that don't pass threshold"""
    if maf_threshold > 0:
        mask_t = maf_t >= maf_threshold
        genotype_t = genotype_t[mask_t]
        variant_id = variant_id[mask_t.cpu().numpy().astype(bool)]
        maf_t = maf_t[mask_t]
    return genotype_t, variant_id,  maf_t