# The script will run tensorqtl with targeted chr (smaller dataset will help speed up the comparison process) 
# Author: Chelsea Trotter 
# Date: Feb 22, 2021

import pandas as pd 
import torch 
import tensorqtl 

from tensorqtl import genotypeio, cis, trans

expression_bed = '../../data/tensorqtldata/GEUVADIS.445_samples.expression.bed.gz'
covariates_file = '../../data/tensorqtldata/GEUVADIS.445_samples.covariates.txt'

# read in genotypes:
chr9_geno_df = pd.read_csv('../../data/tensorqtldata/chr9.csv')  
chr9_geno_df = chr9_geno_df.drop(chr9_geno_df.columns[[1]], axis=1) 
# load phenotype and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

trans_df = trans.map_trans(chr9_geno_df, phenotype_df, covariates_df, batch_size=10000, 
                           return_sparse=True, pval_threshold=1e-5, maf_threshold=0.05)


