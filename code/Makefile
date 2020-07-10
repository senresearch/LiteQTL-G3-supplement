all: getdata installrlib cleandata rscan #compare

.PHONY: all

getdata:
	julia download-data.jl

installrlib:
	Rscript rlibdep.R

cleandata:
	Rscript clean-spleen-data.R
	# Rscript clean-hippo-data.R

# calc_genoprob:
# 	Rscript calc_genoprob.R

# remove_missing_pheno:
# 	Rscript remove_spleen_pheno.R
# 	# Rscript remove_hippo_pheno.R 

rscan:
	# running the scan will also generate geno prob and pheno no missing data file. 
	Rscript Rqtl2scan_spleen.R
	# Rscript Rqtl2scan_hippo.R

#TODO:
# juliascan:
	run_lmgpu.sh
# compare:
	julia check_results.jl
# gemm:
	julia gemm_test.jl
