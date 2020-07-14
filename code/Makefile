all: getdata installlib cleandata rscan #compare

.PHONY: all

getdata:
	julia download-data.jl

getgmap:
	julia getgmap.jl

installlib:
	Rscript rlibdep.R
	julia juliadep.jl

cleandata:
	Rscript clean-spleen-data.R
	# Rscript clean-hippo-data.R

calc_genoprob:
	Rscript calcgenoprob.R

remove_missing_pheno:
	Rscript remove_spleen_pheno.R
# 	# Rscript remove_hippo_pheno.R 

rscan:
	# running the scan will also generate geno prob and pheno no missing data file. 
	Rscript Rqtl2scan_spleen.R
	# Rscript Rqtl2scan_hippo.R

# TODO:
juliascan:
	julia julia-genomescan.jl
# compare:
	julia check_results.jl
# gemm: # TODO: wait until have access to acf. 
# 	julia gemm_test.jl
