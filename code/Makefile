all: rscan juliascan compare
.PHONY: all

getdata:
	julia download-data.jl

installlib:
	# Rscript rlibdep.R
	julia juliadep.jl

#TODO: install CUDA 

cleandata:
	Rscript clean-spleen-data.R
	Rscript clean-hippo-data.R # TODO: running this job on Bluff. Check on this later. 

############## NO LONGER NEEDED ###############
# getgmap:
# 	julia getgmap.jl

# calc_genoprob:
# 	Rscript calcgenoprob.R

# remove_missing_pheno:
# 	Rscript remove_spleen_pheno.R
# 	# Rscript remove_hippo_pheno.R 
##############################################

rscan:
	# running the scan will also generate geno prob and pheno no missing data file. 
	# Rscript Rqtl2scan_spleen.R
	Rscript Rqtl2scan_hippo.R

juliascan:
	JULIA_NUM_THREADS=16 julia julia-genomescan.jl 
	 # runs abount 2 seconds  for max lod and float64 on my home imac 
compare:
	julia check_results.jl
gemm: 
	JULIA_NUM_THREADS=16 julia dgemm-test.jl
	# TODO: generate figure
	julia generate_figure.jl

