import tensorqtl
from tensorqtl import genotypeio

plink_prefix_path = '../../data/tensorqtldata/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup'

pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
# check the first of dataframe
genotype_df.head()
# check the column names of dataframe
list(genotype_df.columns)

# get the row names :
rownames = genotype_df.index.values

# spit chr number from the rest of the row names 
# deliminator is '_', and max split is 1
#chrrownames = [i.split('_', 1) for i in rownames] 
# change data type from object to string 
#strrownames = rownames.astype('str') 
# split and then 
# genotype_df['chr'] = rownames.astype('str') 
# split the string into chr number and the rest of the string 
chr = np.char.split(rownames.astype('str'), sep ='_', maxsplit = 1)
# select only the first element after spliting string 
genotype_df['chr'] = [ i[0] for i in chr ]
# select all the genotype with chr x 
genotype_df.loc[genotype_df['chr'] == 'chrX'] 




