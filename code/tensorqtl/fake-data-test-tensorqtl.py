import numpy as np
import torch 
from scipy import stats
import statistics 

from functionsintensorqtl import *

# r = np.array([[0.1, 0.2, 0.3],  [0.4, 0.5, 0.6]])
# r = np.array([[0.1, 0.2, 0.3],  [0.4, 0.5, 0.6]])

# geno = [0.1 0.8 0.2; 0.6 0.1 0.3]

# pheno = [1.9 1.2 9.3 2.4; 8.5 3.6 9.7 2.8]

geno = np.array([[0.1, 0.8, 10.0, 3.3], [0.2, 0.6, 5.6, 9.3], [2.1, 0.3, 4.1, 9.4], [1.5, 2.4, 8.5, 3.8]])
geno = torch.from_numpy(geno)
pheno = np.array([[1.9, 8.5, 0.6, 0.1], [1.2, 3.6, 2.4, 6.7], [9.3, 9.7, 3.3, 5.8], [2.4, 2.8, 4.5, 9.2]])
pheno = torch.from_numpy(pheno)

dof = geno.size(0) -1

r = calculate_corr(geno, pheno)

res = pval_calc(r, dof)
print(res)