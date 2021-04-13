# The script will run tensorqtl with targeted chr (smaller dataset will help speed up the comparison process) 
# Author: Chelsea Trotter 
# Date: Feb 22, 2021

import pandas as pd 
import torch 
import numpy
import tensorqtl 
import timeit 
import time 
import datetime


from tensorqtl import genotypeio, cis, trans
# from functionsintensorqtl import *

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

small_geno_df = chr9_geno_df.iloc[0:20000, :]
# small_geno_df = chr9_geno_df

##################################### start profiling ########################################

##################################### Full Matrix Case ########################################
for numthreads in [20, 40, 80]:

    torch.set_num_threads(numthreads)
    device = torch.device("cpu")
    trans_df, cpucalctime = trans.map_trans(small_geno_df, phenotype_df, batch_size=20000, 
                            return_sparse=False, pval_threshold=1, maf_threshold=0.00, device=device)

    device = torch.device("cuda")
    (trans_df, gpucalctime) = trans.map_trans(small_geno_df, phenotype_df, batch_size=20000, 
                            return_sparse=False, pval_threshold=1, maf_threshold=0.00, device=device)

    n = small_geno_df.shape[1]
    m = phenotype_df.shape[0]
    p = small_geno_df.shape[0]
    with open("/home/xiaoqihu/git/LiteQTL-G3-supplement/code/tensorqtl/tensorqtl_timing_report.txt", 'a') as file1:
        file1.write(f'{datetime.datetime.now()} \t n:{n},m:{m},p:{p} TensorQTL: \t nthreads: {torch.get_num_threads()} \t cpu(s):{cpucalctime} \t gpu(s){gpucalctime}\n\n')

# trans_df.to_csv('tensorqtl-scan-pval.csv')

##################################### Sparse Matrix Case ########################################
for numthreads in [20,40,80]:

    torch.set_num_threads(numthreads)
    device = torch.device("cpu")
    trans_df, cpucalctime = trans.map_trans(small_geno_df, phenotype_df, batch_size=20000, 
                            return_sparse=True, pval_threshold=1e-5, maf_threshold=0.05, device=device)

    device = torch.device("cuda")
    (trans_df, gpucalctime) = trans.map_trans(small_geno_df, phenotype_df, batch_size=20000, 
                            return_sparse=True, pval_threshold=1e-5, maf_threshold=0.05, device=device)

    n = small_geno_df.shape[1]
    m = phenotype_df.shape[0]
    p = small_geno_df.shape[0]
    with open("/home/xiaoqihu/git/LiteQTL-G3-supplement/code/tensorqtl/tensorqtl_timing_report.txt", 'a') as file1:
        file1.write(f'{datetime.datetime.now()} \t n:{n},m:{m},p:{p} TensorQTL: \t nthreads: {torch.get_num_threads()} \t cpu(s):{cpucalctime} \t gpu(s){gpucalctime}\n\n')



# start_time =  timeit.default_timer()
# trans_df = trans.map_trans(chr9_geno_df, phenotype_df, covariates_df, batch_size=20000, 
#                            return_sparse=True, pval_threshold=1e-5, maf_threshold=0.05)
# #print("Tensorqtl trans.map_trans function took:"
# #print(timeit.default_timer() - start_time)
# timetaken = timeit.default_timer() - start_time
# msg = "{func} took {time} seconds to complete."
# print(msg.format(func = trans.map_trans.__name__, time = timetaken))
# # map_trans took 4.8016955852508545 seconds to complete. # use GPU Nvidia V100

# trans_df.to_csv('tensorqtl-scan-covar.csv')

# import sys
# sys.getsizeof(trans_df) / (1024 *1024 *1024)
