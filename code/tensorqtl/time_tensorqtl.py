# The script will run tensorqtl with targeted chr (smaller dataset will help speed up the comparison process) 
# Author: Chelsea Trotter 
# Date: Feb 22, 2021

import pandas as pd 
import torch 
import numpy
import tensorqtl 

from tensorqtl import genotypeio, cis, trans

expression_bed = '../../data/tensorqtldata/GEUVADIS.445_samples.expression.bed.gz'
covariates_file = '../../data/tensorqtldata/GEUVADIS.445_samples.covariates.txt'

# read in genotypes:
chr9_geno_df = pd.read_csv('../../data/tensorqtldata/chr9.csv')  
chr9_geno_df = chr9_geno_df.drop(chr9_geno_df.columns[[1]], axis=1)
chr9_geno_df = chr9_geno_df.set_index('ID')
#chr9_geno_df = torch.from_numpy(chr9_geno_df.values)
# load phenotype and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

import timeit 
import time 
start_time =  timeit.default_timer()
trans_df = trans.map_trans(chr9_geno_df, phenotype_df, covariates_df, batch_size=10000, 
                           return_sparse=True, pval_threshold=1e-5, maf_threshold=0.05)
#print("Tensorqtl trans.map_trans function took:"
#print(timeit.default_timer() - start_time)
timetaken = timeit.default_timer() - start_time
msg = "{func} took {time} seconds to complete."
print(msg.format(func = trans.map_trans.__name__, time = timetaken))

# trans_df.to_csv('tensorqtl-scan-covar.csv')
# map_trans took 4.8016955852508545 seconds to complete. # use GPU Nvidia V100

# subset data set, getting a smaller end result for the purpose of comparing TensorQTL vs LiteQTL. 
