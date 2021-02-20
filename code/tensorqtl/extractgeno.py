import tensorqtl
from tensorqtl import genotypeio

plink_prefix_path = '../../data/tensorqtldata/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup'

pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()

